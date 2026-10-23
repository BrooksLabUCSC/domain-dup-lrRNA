"""Exon-level tandem duplications and inversions in long-read RNA, in one sample or a tumor-normal pair.

A tandem duplication of exons i..j gives reads whose exons run ... i-1, i .. j, i .. j, j+1 ...: a step back
from exon j to exon i, and a triplication steps back twice, ... i .. j, i .. j, i .. j, j+1 ... An inversion
gives a run of exons in reverse complement inside a read whose other exons are forward. Genome aligners show
these reads in several ways (a large insertion, a supplementary alignment in the same gene, or only a long
soft clip), so the reads are not read from their alignments: each read with one of these signatures has its exons found in its own sequence (edlib), giving the read's exon order.

For each gene (its MANE Select transcript, else its Ensembl canonical one, else its longest coding one):
  1. Reads: primary alignments in the gene on its strand (the reads must be oriented, as Iso-Seq and Kinnex
     reads are). A read has a signature if it has an insertion of at least --min_insertion bases, a soft
     clip of at least --min_clip, or a supplementary alignment in the same gene.
  2. Exon order: each exon of at least --min_exon bases is searched for in the read's sequence, forward and
     reverse complement, allowing --max_error edits per base, every non-overlapping copy.
  3. Events, per read:
       duplication of exons i..j: a forward exon j followed by forward exon i (i <= j); a read stepping back
         from j to i n times shows n + 1 tandem copies of exons i..j
       inversion of exons i..j: a run of reverse complement exons i..j, in a read with forward exons too
     A read supports the event if it also has a forward exon outside i..j. Reads whose exons all lie inside
     i..j (a circular RNA's back-splice gives the same step back) are counted apart as unflanked.
     A duplication is reported once per number of copies, 2 (a duplication), 3 (a triplication) and so on,
     each counting the reads showing at least that many copies: a read of three copies supports both the
     copies 2 and the copies 3 rows, so the copies 2 row counts every read of the event. A read stopping
     before the next copy shows fewer copies than its molecule has.
  4. Reference reads: reads without the event whose alignment has the intron after exon j (or before exon i
     when j is the last exon). fraction = supporting / (supporting + reference).
  5. The event's effect on the protein: a duplication inside the coding sequence is in frame when its length
     is a multiple of 3, with the duplicated residues given; the transcript with its copies is translated from
     the annotated start codon and written to PREFIX.proteins.fa with the reference protein.

An event is reported when at least --min_reads reads of --bam support it. With --normal, each event is also
counted in the normal, which is not searched for events of its own. With --events (an earlier events.tsv),
nothing is searched: those events are counted in --bam, so they can be counted in any other samples.

Outputs:
  PREFIX.events.tsv    one row per event: gene, transcript, event, exons, copies (of a duplication; 1 for an
                       inversion), genomic span, length, coding and frame, duplicated residues, and per sample
                       (sample_, normal_) the supporting, unflanked and reference reads and the fraction
  PREFIX.reads.tsv     each read with an event: its exon order (exon numbers, r for reverse complement)
  PREFIX.proteins.fa   the reference protein and each coding duplication's protein

usage: find_domain_dups.py --bam BAM [--normal BAM] --genome FA --gtf GTF [--genes FILE] -o PREFIX
       find_domain_dups.py --bam BAM --events EVENTS.tsv --genome FA --gtf GTF -o PREFIX
Needs pysam and edlib.
"""
import argparse
import collections
import csv
import multiprocessing
import os
import re

import edlib
import pysam

COMP = str.maketrans('ACGTNacgtn', 'TGCANtgcan')
CODONS = {a + b + c: aa for (a, b, c), aa in zip(
    ((x, y, z) for x in 'TCAG' for y in 'TCAG' for z in 'TCAG'),
    'FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG')}
INTRON_SLOP = 5
COUNT_FIELDS = ['reads', 'unflanked_reads', 'reference_reads', 'fraction']


def revcomp(seq):
    return seq.translate(COMP)[::-1]


def translate(seq):
    "protein up to the first stop codon"
    out = []
    for i in range(0, len(seq) - 2, 3):
        aa = CODONS.get(seq[i:i + 3].upper(), 'X')
        if aa == '*':
            break
        out.append(aa)
    return ''.join(out)


class Gene:
    """a gene's chosen transcript: exons (genomic, in transcript order), their sequences and the CDS"""

    def __init__(self, gene_id, name, chrom, strand, start, end, transcript, exons, cds):
        self.id, self.name, self.chrom, self.strand = gene_id, name, chrom, strand
        self.start, self.end, self.transcript = start, end, transcript
        self.exons = sorted(exons, reverse=strand == '-')
        self.cds = cds   # genomic (start, end) of the CDS, None if non-coding
        self.seqs = None

    def load(self, genome):
        seqs = [genome.fetch(self.chrom, s, e).upper() for s, e in self.exons]
        self.seqs = [revcomp(x) for x in seqs] if self.strand == '-' else seqs
        offsets = [0]
        for s in self.seqs:
            offsets.append(offsets[-1] + len(s))
        self.offsets = offsets   # transcript coordinate of each exon's first base, and the length last
        self.cds_tx = None
        if self.cds:
            first = self.cds[0] if self.strand == '+' else self.cds[1] - 1
            last = self.cds[1] - 1 if self.strand == '+' else self.cds[0]
            self.cds_tx = (self.tx_pos(first), self.tx_pos(last) + 1)

    def tx_pos(self, pos):
        for k, (s, e) in enumerate(self.exons):
            if s <= pos < e:
                return self.offsets[k] + (pos - s if self.strand == '+' else e - 1 - pos)
        raise ValueError(f'{pos} not in an exon of {self.transcript}')

    def intron(self, k):
        "genomic (start, end) of the intron after exon k (0-based, transcript order)"
        a, b = self.exons[k], self.exons[k + 1]
        return (a[1], b[0]) if self.strand == '+' else (b[1], a[0])


def attributes(text):
    out, tags = {}, []
    for key, value in re.findall(r'(\S+) "([^"]*)"', text):
        if key == 'tag':
            tags.append(value)
        else:
            out[key] = value
    out['tags'] = tags
    return out


def field(text, key):
    "one attribute's value, without parsing the rest"
    i = text.find(key + ' "')
    if i < 0:
        return None
    i += len(key) + 2
    return text[i:text.index('"', i)]


def read_genes(gtf, wanted):
    """{gene_id: Gene} for the wanted genes (names or ids, versions optional), or every protein-coding gene"""
    genes, transcripts = {}, collections.defaultdict(lambda: {'exons': [], 'cds': [], 'tags': []})
    for line in open(gtf):
        if line.startswith('#'):
            continue
        f = line.rstrip('\n').split('\t')
        if f[2] not in ('gene', 'transcript', 'exon', 'CDS'):
            continue
        gid = field(f[8], 'gene_id')
        name = field(f[8], 'gene_name') or gid
        if wanted is not None:
            if name not in wanted and gid not in wanted and gid.split('.')[0] not in wanted:
                continue
        elif (field(f[8], 'gene_type') or field(f[8], 'gene_biotype') or 'protein_coding') != 'protein_coding':
            continue
        a = attributes(f[8])
        if f[2] == 'gene':
            genes[gid] = (name, f[0], f[6], int(f[3]) - 1, int(f[4]))
            continue
        t = transcripts[(gid, a['transcript_id'])]
        if f[2] == 'transcript':
            t['tags'] = a['tags']
        else:
            t['exons' if f[2] == 'exon' else 'cds'].append((int(f[3]) - 1, int(f[4])))
    by_gene = collections.defaultdict(list)
    for (gid, tid), t in transcripts.items():
        by_gene[gid].append((tid, t))
    out = {}
    for gid, ts in by_gene.items():
        if gid not in genes or not ts:
            continue

        def rank(item):
            tid, t = item
            return ('MANE_Select' in t['tags'], 'Ensembl_canonical' in t['tags'],
                    sum(e - s for s, e in t['cds']), sum(e - s for s, e in t['exons']), tid)
        tid, t = max(ts, key=rank)
        name, chrom, strand, start, end = genes[gid]
        cds = (min(s for s, _ in t['cds']), max(e for _, e in t['cds'])) if t['cds'] else None
        out[gid] = Gene(gid, name, chrom, strand, start, end, tid, t['exons'], cds)
    return out


def signature(a, gene, args):
    "whether an alignment shows a rearrangement: a large insertion or soft clip, or a split in the gene"
    for op, n in a.cigartuples:
        if (op == 1 and n >= args.min_insertion) or (op == 4 and n >= args.min_clip):
            return True
    if a.has_tag('SA'):
        for sa in a.get_tag('SA').rstrip(';').split(';'):
            chrom, pos = sa.split(',')[:2]
            if chrom == gene.chrom and gene.start <= int(pos) - 1 < gene.end:
                return True
    return False


def exon_hits(seq, gene, args):
    """[(read start, read end, exon number (1-based), strand)] of the exons found in a read, in read order,
    overlaps resolved in favor of the fewer edits per base"""
    hits = []
    for k, ex in enumerate(gene.seqs):
        if len(ex) < args.min_exon:
            continue
        limit = max(1, int(args.max_error * len(ex)))
        for strand, query in (('+', ex), ('-', revcomp(ex))):
            masked = seq
            for _ in range(args.max_copies):
                r = edlib.align(query, masked, mode='HW', task='locations', k=limit)
                if r['editDistance'] < 0:
                    break
                taken = []
                for s, e in r['locations']:
                    s = s or 0
                    if all(e < ts or s > te for ts, te in taken):
                        taken.append((s, e))
                for s, e in taken:
                    hits.append((s, e + 1, k + 1, strand, r['editDistance'] / len(ex)))
                    masked = masked[:s] + 'N' * (e + 1 - s) + masked[e + 1:]
    hits.sort(key=lambda h: (h[4], h[0]))
    kept = []
    for h in hits:   # best first: drop a hit that overlaps a kept one by more than half the shorter
        if all(min(h[1], k[1]) - max(h[0], k[0]) <= 0.5 * min(h[1] - h[0], k[1] - k[0]) for k in kept):
            kept.append(h)
    return sorted((h[:4] for h in kept), key=lambda h: h[0])


def read_events(hits):
    """{(event, i, j, copies): flanked} of a read's exon order.  A duplication of exons i..j the read steps
    back over n times shows n + 1 copies, and supports the events of 2 to n + 1 copies; an inversion has 1"""
    fwd = [h[2] for h in hits if h[3] == '+']
    steps = collections.Counter((b, a) for a, b in zip(fwd, fwd[1:]) if b <= a)
    out = {}
    for (i, j), n in steps.items():
        flanked = any(x < i or x > j for x in fwd)
        for copies in range(2, n + 2):
            out[('duplication', i, j, copies)] = flanked
    if fwd:
        run = []
        for h in hits + [(0, 0, 0, '+')]:
            if h[3] == '-':
                run.append(h[2])
            elif run:
                i, j = min(run), max(run)
                out[('inversion', i, j, 1)] = any(x < i or x > j for x in fwd)
                run = []
    return out


def has_intron(a, intron):
    pos = a.reference_start
    for op, n in a.cigartuples:
        if op == 3 and abs(pos - intron[0]) <= INTRON_SLOP and abs(pos + n - intron[1]) <= INTRON_SLOP:
            return True
        if op in (0, 2, 3, 7, 8):
            pos += n
    return False


def scan(bam, gene, args):
    """the reads of a gene: ({event: [supporting, unflanked]}, {read: exon order}, alignments without events)"""
    counts = collections.defaultdict(lambda: [0, 0])
    orders, plain = {}, []
    with pysam.AlignmentFile(bam) as f:
        for a in f.fetch(gene.chrom, gene.start, gene.end):
            if a.is_secondary or a.is_supplementary or a.is_unmapped or a.is_reverse != (gene.strand == '-'):
                continue
            ev = {}
            if signature(a, gene, args):
                seq = a.get_forward_sequence()
                if seq:
                    hits = exon_hits(seq.upper(), gene, args)
                    ev = read_events(hits)
                    if ev:
                        orders[a.query_name] = ' '.join(f"{h[2]}{'r' if h[3] == '-' else ''}" for h in hits)
            for key, flanked in ev.items():
                counts[key][0 if flanked else 1] += 1
            if not any(ev.values()):
                plain.append(a)
    return counts, orders, plain


def reference_reads(plain, gene, event):
    _, i, j, _ = event
    k = j - 1 if j < len(gene.exons) else i - 2
    if k < 0 or k + 1 >= len(gene.exons):
        return None
    intron = gene.intron(k)
    return sum(has_intron(a, intron) for a in plain)


def effect(gene, event):
    """(coding, in_frame, residues, protein) of a duplication of exons i..j, in its number of copies"""
    kind, i, j, copies = event
    if kind != 'duplication' or gene.cds_tx is None:
        return 'no', '', '', None
    s, e = gene.offsets[i - 1], gene.offsets[j]
    cs, ce = gene.cds_tx
    if e <= cs or s >= ce:
        return 'no', '', '', None
    tx = ''.join(gene.seqs)
    model = tx[:e] + tx[s:e] * (copies - 1) + tx[e:]
    protein = translate(model[cs:])
    if s < cs or e > ce:
        return 'partly', 'no', '', protein
    in_frame = (e - s) % 3 == 0
    residues = f'{(s - cs) // 3 + 1}-{(e - cs - 1) // 3 + 1}' if in_frame else ''
    return 'yes', 'yes' if in_frame else 'no', residues, protein


def span(gene, i, j):
    blocks = gene.exons[i - 1:j]
    return min(s for s, _ in blocks), max(e for _, e in blocks)


def run_gene(job):
    gene, args, wanted = job
    gene.load(pysam.FastaFile(args.genome))
    counts, orders, plain = scan(args.bam, gene, args)
    if wanted is None:
        events = [k for k, v in counts.items() if v[0] >= args.min_reads]
    else:
        events = wanted
    if not events:
        return [], [], []
    samples = [('sample', counts, plain)]
    if args.normal:
        n_counts, _, n_plain = scan(args.normal, gene, args)
        samples.append(('normal', n_counts, n_plain))
    rows, proteins = [], []
    if gene.cds_tx:
        tx = ''.join(gene.seqs)
        proteins.append((f'{gene.name}|{gene.transcript}|reference', translate(tx[gene.cds_tx[0]:])))
    for ev in sorted(events, key=lambda x: (x[1], x[2], x[0], x[3])):
        kind, i, j, copies = ev
        start, end = span(gene, i, j)
        coding, in_frame, residues, protein = effect(gene, ev)
        exons = f'{i}-{j}' if i != j else str(i)
        row = {'gene': gene.name, 'gene_id': gene.id, 'transcript': gene.transcript, 'event': kind, 'exons': exons,
               'copies': copies, 'exon_count': len(gene.exons), 'chrom': gene.chrom, 'start': start, 'end': end,
               'strand': gene.strand,
               'length': gene.offsets[j] - gene.offsets[i - 1], 'coding': coding, 'in_frame': in_frame,
               'residues': residues, 'protein_length': len(protein) if protein else ''}
        for label, c, p in samples:
            sup, unf = c.get(ev, (0, 0))
            ref = reference_reads(p, gene, ev)
            row.update({f'{label}_reads': sup, f'{label}_unflanked_reads': unf,
                        f'{label}_reference_reads': '' if ref is None else ref,
                        f'{label}_fraction': '' if not ref and not sup else f'{sup / (sup + (ref or 0)):.4f}'})
        rows.append(row)
        if protein:
            proteins.append((f'{gene.name}|{gene.transcript}|{kind}_{exons}' + (f'_x{copies}' if copies > 2 else ''),
                             protein))
    read_rows = [(r, gene.name, o) for r, o in orders.items()]
    return rows, read_rows, proteins


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--bam', required=True, help='the sample (the tumor of a pair), indexed')
    parser.add_argument('--normal', help='a matched normal, whose reads are counted for every event')
    parser.add_argument('--events', help='count these events (an events.tsv) instead of searching')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--gtf', required=True)
    parser.add_argument('--genes', help='genes to search, a name or id per line, or a BED with names in column 4 '
                        '(default every protein-coding gene)')
    parser.add_argument('-o', '--prefix', required=True)
    parser.add_argument('--sample_name', help='written in the sample_name column (default the BAM name)')
    parser.add_argument('--normal_name')
    parser.add_argument('--min_reads', type=int, default=3, help='supporting reads in --bam to report an event')
    parser.add_argument('--min_insertion', type=int, default=30)
    parser.add_argument('--min_clip', type=int, default=30)
    parser.add_argument('--min_exon', type=int, default=12, help='shorter exons are not searched for')
    parser.add_argument('--max_error', type=float, default=0.05, help='edits per base allowed in an exon match')
    parser.add_argument('--max_copies', type=int, default=6, help='copies of an exon searched for in a read')
    parser.add_argument('-t', '--threads', type=int, default=1)
    args = parser.parse_args()

    wanted_events = None
    if args.events:
        wanted_events = collections.defaultdict(list)
        for r in csv.DictReader(open(args.events), delimiter='\t'):
            i, _, j = r['exons'].partition('-')
            # an events.tsv from before copies were counted is of duplications of 2 copies
            copies = int(r.get('copies') or (2 if r['event'] == 'duplication' else 1))
            wanted_events[r['gene_id']].append((r['event'], int(i), int(j or i), copies))
        names = set(wanted_events)
    elif args.genes:
        names = {line.rstrip('\n').split('\t')[3] if line.count('\t') >= 3 else line.strip()
                 for line in open(args.genes) if line.strip() and not line.startswith('#')}
    else:
        names = None
    genes = read_genes(args.gtf, names)
    if wanted_events is not None:
        jobs = [(g, args, sorted(set(wanted_events[g.id]))) for g in genes.values()]
    else:
        jobs = [(g, args, None) for g in genes.values()]
    jobs.sort(key=lambda x: (x[0].chrom, x[0].start))

    sample = args.sample_name or os.path.basename(args.bam).split('.')[0]
    normal = (args.normal_name or os.path.basename(args.normal).split('.')[0]) if args.normal else ''
    rows, reads, proteins = [], [], []
    with multiprocessing.Pool(args.threads) as pool:
        for r, rd, p in pool.imap(run_gene, jobs, chunksize=4):
            rows.extend(r)
            reads.extend(rd)
            proteins.extend(p)

    labels = ['sample'] + (['normal'] if args.normal else [])
    fields = (['sample_name', 'normal_name', 'gene', 'gene_id', 'transcript', 'event', 'exons', 'copies', 'exon_count',
               'chrom',
               'start', 'end', 'strand', 'length', 'coding', 'in_frame', 'residues', 'protein_length']
              + [f'{lab}_{x}' for lab in labels for x in COUNT_FIELDS])
    with open(f'{args.prefix}.events.tsv', 'w') as f:
        w = csv.DictWriter(f, fields, delimiter='\t', lineterminator='\n', restval='')
        w.writeheader()
        for r in rows:
            w.writerow(dict(r, sample_name=sample, normal_name=normal))
    with open(f'{args.prefix}.reads.tsv', 'w') as f:
        f.write('sample_name\tread\tgene\texon_order\n')
        for r, g, o in reads:
            f.write(f'{sample}\t{r}\t{g}\t{o}\n')
    with open(f'{args.prefix}.proteins.fa', 'w') as f:
        for name, seq in proteins:
            f.write(f'>{name}\n{seq}\n')
    print(f'{len(genes)} genes, {len(rows)} events. Wrote {args.prefix}.events.tsv')


if __name__ == '__main__':
    main()
