"""Transcript models of the events find_domain_dups.py found in a gene, and every read over them realigned.

find_domain_dups.py finds duplications and inversions from each read's exon order, without aligning reads
across them. To call and phase variants on a duplicated region, as on the copies of a kinase domain, the reads
need alignments to the transcripts with the event. For one gene of an events.tsv, this builds:

  models   the gene's transcript (the one find_domain_dups.py used, from the same GTF), exons 1..N, as the
           normal model, and a model of each of the gene's events in the events.tsv: a duplication of exons
           i..j in c copies is 1..j, then i..j c - 1 more times, then j+1..N (copies 2 a duplication, 3 a
           triplication), and an inversion is 1..i-1, j..i reverse complemented, j+1..N
  reads    every read with an alignment over any of the events' exons, in each --bam, not only the reads
           find_domain_dups.py found the events in, and including reads whose primary alignment is elsewhere
           in the gene. Each read's primary alignment gives its sequence and base qualities, 5' to 3', and
           its --keep_tags (default HP and PS, the haplotags of whatshap haplotag), which minimap2 -y copies
           to the new alignment
  realign  the reads aligned to the models with minimap2 (--secondary=no --MD -y), sorted and indexed

A read's model is minimap2's choice. A read that doesn't reach an event's junctions aligns equally well to
several models and gets a mapping quality of 0, so each event is counted by the junctions reads cross, as
find_domain_dups.py counts it from exon orders, whatever the model and mapping quality:
  duplication of exons i..j in c copies: reads crossing the junction from exon j to exon i at least c - 1
    times, so showing at least c copies (a read of three copies counts for the copies 2 and 3 events)
  inversion of exons i..j: reads crossing into or out of the reverse complemented exons
  normal: reads of the normal model, of at least --min_mapq, crossing the junction the event replaces (exon
    j to j+1, or i-1 to i when j is the last exon)
fraction = event reads / (event reads + normal reads), as find_domain_dups.py takes reference reads.

Variants can then be called on PREFIX.<sample>.bam with PREFIX.models.fa as the reference, and phased with
the haplotags the reads carry.

Outputs:
  PREFIX.models.fa (.fai)     the models, named GENE_normal, GENE_duplication_18-25, GENE_duplication_18-25_x3,
                              GENE_inversion_5-6 and so on
  PREFIX.models.exons.bed     each model's exons in model coordinates, named by exon number, with .2, .3 for
                              further copies and r for reverse complement (load with the BAMs in IGV)
  PREFIX.<sample>.bam (.bai)  each sample's reads aligned to the models
  PREFIX.events.tsv           per sample and event: its model, event and normal reads, and fraction
  PREFIX.models.tsv           per sample and model: reads, and reads of at least --min_mapq
  PREFIX.reads.tsv            per read: sample, model, mapping quality, the junctions it crosses, the events
                              it supports

usage: build_dup_models.py --events T1.events.tsv --gene EGFR --bam T1.bam [--bam N1.bam] --genome FA --gtf GTF
           -o PREFIX
Needs pysam, and minimap2 and samtools in PATH. Imports find_domain_dups.py from its own directory.
"""
import argparse
import collections
import csv
import os
import subprocess
import sys
import tempfile

import pysam

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from find_domain_dups import read_genes, revcomp  # noqa: E402

JUNCTION_FLANK = 20

Event = collections.namedtuple('Event', ('kind', 'i', 'j', 'copies', 'chrom', 'start', 'end'))


def gene_events(events_tsv, gene):
    """the gene's events of an events.tsv, each once, and its gene id and transcript"""
    events, ids = {}, set()
    for r in csv.DictReader(open(events_tsv), delimiter='\t'):
        if gene not in (r['gene'], r['gene_id'], r['gene_id'].split('.')[0]):
            continue
        ids.add((r['gene_id'], r['transcript']))
        i, _, j = r['exons'].partition('-')
        copies = int(r.get('copies') or (2 if r['event'] == 'duplication' else 1))
        key = (r['event'], int(i), int(j or i), copies)
        events[key] = Event(*key, r['chrom'], int(r['start']), int(r['end']))
    if not events:
        sys.exit(f'{events_tsv} has no events of {gene}')
    if len(ids) > 1:
        sys.exit(f'{events_tsv}: {gene} is of more than one gene or transcript: {sorted(ids)}')
    return sorted(events.values(), key=lambda e: (e.i, e.j, e.kind, e.copies)), ids.pop()


def model_name(gene, ev):
    return f'{gene.name}_{ev.kind}_{ev.i}-{ev.j}' + (f'_x{ev.copies}' if ev.copies > 2 else '')


def exon_order(ev, n):
    "an event's model exon order, [(exon, strand)]"
    if ev.kind == 'duplication':
        return ([(k, '+') for k in range(1, ev.j + 1)]
                + [(k, '+') for _ in range(ev.copies - 2) for k in range(ev.i, ev.j + 1)]
                + [(k, '+') for k in range(ev.i, n + 1)])
    return ([(k, '+') for k in range(1, ev.i)] + [(k, '-') for k in range(ev.j, ev.i - 1, -1)]
            + [(k, '+') for k in range(ev.j + 1, n + 1)])


def build_models(gene, events):
    """[(model name, sequence, [(start, end, label, exon, strand)])]: the normal model, then one per event"""
    n = len(gene.exons)
    named = [(f'{gene.name}_normal', [(k, '+') for k in range(1, n + 1)])]
    named += [(model_name(gene, ev), exon_order(ev, n)) for ev in events]
    models = []
    for name, order in named:
        seq, exons, seen, pos = [], [], collections.Counter(), 0
        for k, strand in order:
            ex = gene.seqs[k - 1] if strand == '+' else revcomp(gene.seqs[k - 1])
            seen[(k, strand)] += 1
            label = f'{k}{"r" if strand == "-" else ""}' + ('' if seen[(k, strand)] == 1 else f'.{seen[(k, strand)]}')
            exons.append((pos, pos + len(ex), label, k, strand))
            seq.append(ex)
            pos += len(ex)
        models.append((name, ''.join(seq), exons))
    return models


def read_names(bam, events):
    "the reads with any alignment (primary or supplementary) over any event's exons"
    names = set()
    with pysam.AlignmentFile(bam) as f:
        for ev in events:
            for a in f.fetch(ev.chrom, ev.start, ev.end):
                if not a.is_secondary and not a.is_unmapped:
                    names.add(a.query_name)
    return names


def write_fastq(bam, gene, names, keep_tags, out):
    """each named read's primary alignment in the gene's region, 5' to 3', with its kept tags as a FASTQ
    comment; returns the reads written"""
    written = set()
    with pysam.AlignmentFile(bam) as f, open(out, 'w') as fq:
        for a in f.fetch(gene.chrom, gene.start, gene.end):
            if a.query_name not in names or a.query_name in written or a.is_secondary or a.is_supplementary:
                continue
            seq = a.get_forward_sequence()
            if not seq:
                continue
            # a read without base qualities, as from FASTA, is given a placeholder of 40 (I)
            quals = a.get_forward_qualities() if a.query_qualities is not None else None
            qual = ''.join(chr(q + 33) for q in quals) if quals is not None else 'I' * len(seq)
            tags = []
            for tag in keep_tags:
                if a.has_tag(tag):
                    value, vtype = a.get_tag(tag, with_value_type=True)
                    tags.append(f'{tag}:{"i" if vtype in "cCsSiI" else vtype}:{value}')
            comment = ('\t' + '\t'.join(tags)) if tags else ''
            fq.write(f'@{a.query_name}{comment}\n{seq}\n+\n{qual}\n')
            written.add(a.query_name)
    return written


def realign(models_fa, fastq, out_bam, args):
    mm2 = ['minimap2', '-a', '--secondary=no', '--MD', '-y', '-t', str(args.threads)]
    if args.preset:
        mm2 += ['-x', args.preset]
    mm2 += [models_fa, fastq]
    p1 = subprocess.Popen(mm2, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL)
    subprocess.run(['samtools', 'sort', '-o', out_bam, '-'], stdin=p1.stdout, check=True)
    p1.stdout.close()
    if p1.wait():
        sys.exit(f'minimap2 failed on {fastq}')
    subprocess.run(['samtools', 'index', out_bam], check=True)


def crossed_junctions(a, exons):
    """the model's exon junctions the alignment crosses, by JUNCTION_FLANK bases each side, as
    [((exon, strand), (next exon, strand))]"""
    out = []
    for x, y in zip(exons, exons[1:]):
        if a.reference_start <= x[1] - JUNCTION_FLANK and a.reference_end >= x[1] + JUNCTION_FLANK:
            out.append(((x[3], x[4]), (y[3], y[4])))
    return out


def supports(ev, crossed):
    "whether a read crossing these junctions shows the event"
    if ev.kind == 'duplication':
        return sum(1 for a, b in crossed if a == (ev.j, '+') and b == (ev.i, '+')) >= ev.copies - 1
    inside = lambda x: x[1] == '-' and ev.i <= x[0] <= ev.j  # noqa: E731
    return any(inside(a) != inside(b) for a, b in crossed)


def normal_junction(ev, n):
    k = ev.j if ev.j < n else ev.i - 1
    return ((k, '+'), (k + 1, '+'))


def label(junction):
    (a, sa), (b, sb) = junction
    return f'{a}{"r" if sa == "-" else ""}-{b}{"r" if sb == "-" else ""}'


def tally(sample, out_bam, gene, events, models, args, read_rows):
    """per event [event reads, normal reads], and per model [reads, reads >= min_mapq]"""
    exons = {name: ex for name, _, ex in models}
    normal = f'{gene.name}_normal'
    n = len(gene.exons)
    event_counts = {ev: [0, 0] for ev in events}
    model_counts = {name: [0, 0] for name in exons}
    with pysam.AlignmentFile(out_bam) as f:
        for a in f:
            if a.is_unmapped or a.is_secondary or a.is_supplementary:
                continue
            crossed = crossed_junctions(a, exons[a.reference_name])
            model_counts[a.reference_name][0] += 1
            model_counts[a.reference_name][1] += a.mapping_quality >= args.min_mapq
            shown = []
            for ev in events:
                if supports(ev, crossed):
                    event_counts[ev][0] += 1
                    shown.append(model_name(gene, ev))
                elif (a.reference_name == normal and a.mapping_quality >= args.min_mapq
                      and normal_junction(ev, n) in crossed):
                    event_counts[ev][1] += 1
            read_rows.append((sample, a.query_name, a.reference_name, a.mapping_quality,
                              ','.join(label(j) for j in crossed), ','.join(shown)))
    return event_counts, model_counts


def run_sample(sample, bam, gene, events, models, models_fa, args, read_rows):
    "one sample's reads over the events, realigned and counted: its events and models rows"
    names = read_names(bam, events)
    keep_tags = [t for t in args.keep_tags.split(',') if t]
    with tempfile.TemporaryDirectory(dir=os.path.dirname(os.path.abspath(args.prefix))) as tmp:
        fastq = os.path.join(tmp, 'reads.fq')
        n = len(write_fastq(bam, gene, names, keep_tags, fastq))
        out_bam = f'{args.prefix}.{sample}.bam'
        realign(models_fa, fastq, out_bam, args)
    event_counts, model_counts = tally(sample, out_bam, gene, events, models, args, read_rows)
    event_rows = []
    for ev, (sup, ref) in event_counts.items():
        exons = f'{ev.i}-{ev.j}' if ev.i != ev.j else str(ev.i)
        fraction = f'{sup / (sup + ref):.4f}' if sup + ref else ''
        event_rows.append((sample, gene.name, ev.kind, exons, ev.copies, model_name(gene, ev), sup, ref, fraction))
    print(f'{sample}: {n} reads realigned; ' + ', '.join(f'{r[5]} {r[6]} event / {r[7]} normal reads'
                                                         for r in event_rows))
    return event_rows, [(sample, m, *c) for m, c in model_counts.items()]


def write_models(models, prefix):
    "the models' FASTA (indexed) and their exons as a BED; returns the FASTA's path"
    models_fa = f'{prefix}.models.fa'
    with open(models_fa, 'w') as fa, open(f'{prefix}.models.exons.bed', 'w') as bed:
        for name, seq, exons in models:
            fa.write(f'>{name}\n{seq}\n')
            for s, e, lab, _, _ in exons:
                bed.write(f'{name}\t{s}\t{e}\t{lab}\n')
    pysam.faidx(models_fa)
    return models_fa


def write_table(path, header, rows):
    with open(path, 'w') as f:
        f.write('\t'.join(header) + '\n')
        for r in rows:
            f.write('\t'.join(map(str, r)) + '\n')


def main():
    parser = argparse.ArgumentParser(description=__doc__.split('\n')[0],
                                     formatter_class=argparse.RawDescriptionHelpFormatter, epilog=__doc__)
    parser.add_argument('--events', required=True, help='an events.tsv of find_domain_dups.py')
    parser.add_argument('--gene', required=True, help='the gene, by name or id: every one of its events is modelled')
    parser.add_argument('--bam', required=True, action='append', help='a sample\'s genome-aligned BAM, indexed; '
                        'give it once per sample, as the tumor and its normal')
    parser.add_argument('--sample_names', help='the samples\' names, comma separated (default the BAM names)')
    parser.add_argument('--genome', required=True)
    parser.add_argument('--gtf', required=True, help='the GTF find_domain_dups.py was run with')
    parser.add_argument('-o', '--prefix', required=True)
    parser.add_argument('--keep_tags', default='HP,PS', help='read tags carried to the new alignments')
    parser.add_argument('--min_mapq', type=int, default=2, help='mapping quality to count a normal read')
    parser.add_argument('--preset', help='a minimap2 -x preset (default none, as the original analysis)')
    parser.add_argument('-t', '--threads', type=int, default=4)
    args = parser.parse_args()

    events, (gene_id, transcript) = gene_events(args.events, args.gene)
    genes = read_genes(args.gtf, {gene_id})
    if gene_id not in genes or genes[gene_id].transcript != transcript:
        sys.exit(f'{args.gtf} doesn\'t give {gene_id} the transcript of the events, {transcript}: give the GTF '
                 'find_domain_dups.py was run with')
    gene = genes[gene_id]
    gene.load(pysam.FastaFile(args.genome))
    models = build_models(gene, events)
    models_fa = write_models(models, args.prefix)

    names = (args.sample_names.split(',') if args.sample_names
             else [os.path.basename(b).split('.')[0] for b in args.bam])
    if len(names) != len(args.bam):
        sys.exit('--sample_names needs a name per --bam')
    read_rows, event_rows, model_rows = [], [], []
    for sample, bam in zip(names, args.bam):
        ev_rows, m_rows = run_sample(sample, bam, gene, events, models, models_fa, args, read_rows)
        event_rows += ev_rows
        model_rows += m_rows
    write_table(f'{args.prefix}.events.tsv', ['sample_name', 'gene', 'event', 'exons', 'copies', 'model',
                                              'event_reads', 'normal_reads', 'fraction'], event_rows)
    write_table(f'{args.prefix}.models.tsv', ['sample_name', 'model', 'reads', f'reads_mapq{args.min_mapq}'],
                model_rows)
    write_table(f'{args.prefix}.reads.tsv', ['sample_name', 'read', 'model', 'mapq', 'junctions', 'events'],
                read_rows)


if __name__ == '__main__':
    main()
