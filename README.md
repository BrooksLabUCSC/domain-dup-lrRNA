# Exon-level domain duplications from long-read RNA

Scripts that find tandem duplications of exons (such as kinase domain duplications)
in long-read RNA sequencing, from the order of exons within each read.

| File | What it is |
|---|---|
| `original/identify_domain_dup_v2.py` | The script used for the published analysis, kept as it was run |
| `find_domain_dups.py` | A cleaned-up, general version of the search, with a command line |
| `build_dup_models.py` | The models of a gene's events from `find_domain_dups.py`, and every read over them realigned to them, for calling and phasing variants on the duplication |

A tandem duplication of exons i..j makes transcripts whose exons run
... i-1, i .. j, i .. j, j+1 ... Genome aligners show such reads in several ways (a
large insertion, a supplementary alignment in the same gene, or a long soft clip), so
both scripts find each candidate read's exons in its own sequence with
[edlib](https://github.com/Martinsos/edlib), and call duplications from the exon order.

## The original analysis (`original/identify_domain_dup_v2.py`)

This is the script as it was run for the paper, on tumor-normal or tumor-only long-read RNA
sequencing (PacBio). It has the paths, sample naming and gene lists of that dataset
written into it, so it is shared to document the analysis rather than to be rerun
as is. `find_domain_dups.py` (below) is the version to run on new data.

### Preprocessing

Each sample's reads were aligned to the genome (GRCh38) with minimap2 in spliced
mode, writing MD tags, giving `<sample>.genomealigned.withMD.bam`. Then, for each
sample:

1. **Split reads.** The primary alignments with a supplementary alignment (an `SA`
   tag) were written to a chimeric-read BAM, converted to BED
   (`bedtools bamtobed`) and intersected with gene bounds from the annotation
   (`bedtools intersect -wo`), giving `<sample>.temp.genomicchim.geneoverlap.txt`.
   A read with two or more alignments in the same gene is evidence of a
   rearrangement within that gene.
2. **Large insertions.** Indels were called from the genome-aligned reads in the
   regions of the cancer genes, at a minimum variant allele fraction of 5%, with
   FLAIR's `identify_vars.py`. For a tumor with a matched normal, indels the normal
   also had were removed. Insertions of at least 100 bases overlapping a MANE Select
   exon (plus 2 bases) were kept in `<sample>.indels.exonOverlap.largeInsertions.bed`
   (`combine_annotate_variants.py`, `identify_large_insertions`), and these were
   intersected with the gene bounds (`bedtools intersect -wo`), giving
   `<sample>.temp.largeInserts.withGenes.txt`.
3. **Reference files.** The genome FASTA, the GENCODE v49 basic GTF (from which the
   MANE Select exons of each gene are taken), a table of Ensembl gene ids to gene
   names, and the genes searched: an OncoPlex v6 panel gene list
   (`oncoplex_v6_gene_list.csv`) and the OncoKB cancer gene list (`cancerGeneList.tsv`).

### Method

For each tumor, tumor recurrence and cell line sample:

1. **Candidate genes**: genes with at least 10 split reads within the gene, plus the
   reads of the gene's large insertions, where these are more than 1% of the gene's
   reads and the gene is in the cancer gene list.
2. **Exon order**: for each primary read on the gene's strand, each MANE Select exon
   of at least 20 bases, forward and reverse complement, is searched for in the
   read's sequence with edlib (at most 5 edits), every copy found. The exons found,
   in read order, are the read's exon order.
3. **Models**: exon orders seen in at least 3 reads, not contained in a longer one,
   that repeat an exon or have one in reverse complement mark the duplicated exons.
   Each is completed to a full transcript (the missing first and last exons added,
   and skipped exons filled in), and written with the normal MANE transcript as a
   small reference of transcript models.
4. **Counting**: the reads of the gene's locus, of the sample and of its matched
   normal, are aligned to the models with minimap2. Reads with less than 50 bases
   of insertion are assigned the exons they cover.
5. **Call**: a duplication is reported when at least 3 of the sample's reads cover a
   duplicated model and no read of the matched normal does. Samples without a
   matched normal (cell lines, recurrences) need only the 3 reads.

Output: `all.dup.final.<date>.txt`, one row per call: sample, gene, duplicated exons
(MANE Select exon numbers), reads supporting the duplication, reads covering the
duplicated exons, and their ratio.

Notes on the archived script:

- The `bedtools` steps that made the overlap files are commented out in it, as they
  were run once beforehand.
- It has a line that skips the first 23 samples (`if c < 24: continue`), left from
  restarting an interrupted run, and it leaves out five genes by name (AFDN, ATRX,
  PSIP1, CCAR1, TOP1).

## The new version (`find_domain_dups.py`)

The same idea, rewritten as a general command-line tool. The differences from the
original:

- No preprocessing: it reads the genome-aligned BAM directly, and picks the reads
  with a large insertion, a long soft clip or a supplementary alignment in the gene
  itself.
- No realignment to transcript models in the search: duplications and inversions are
  called from each read's exon order, and a read supports an event only if it also
  has an exon outside the duplicated ones, so circular RNA back-splices are counted
  apart. This makes it fast enough to search every gene. The models and realignment
  of the original are a separate step, `build_dup_models.py` (below), run on the
  genes of interest.
- Copies: a read whose exon order steps back from exon j to exon i twice shows three
  tandem copies of exons i..j. A duplication is reported once per number of copies
  with at least `--min_reads` reads (copies 2, a duplication, copies 3, a
  triplication, and so on), each counting the reads showing at least that many
  copies, so the copies 2 row counts every read of the event.
- A fraction: supporting reads against reference reads, those spanning the intron
  after the duplicated exons without the event.
- The effect on the protein: whether a duplication is in the coding sequence and in
  frame, the duplicated residues, and the duplicated protein.
- Any genes (default every protein-coding gene), and any samples: a matched normal is
  counted for every event found, and a table of events can be counted in other
  samples, such as a panel of normals.

### Requirements

Python 3 (tested with 3.11) with [pysam](https://github.com/pysam-developers/pysam) and
[edlib](https://github.com/Martinsos/edlib) (`pip install pysam edlib`, or from
conda-forge and bioconda).

Inputs:

- reads aligned to the genome, as a sorted and indexed BAM, oriented to their RNA's
  strand (as PacBio Iso-Seq and Kinnex reads are). Align with minimap2 in spliced
  mode, for example
  `minimap2 -ax splice -s 40 -G 350k --secondary=no GENOME.fa reads.fa | samtools sort -o sample.bam - && samtools index sample.bam`
- the genome FASTA (indexed with `samtools faidx`) and a GTF annotation. A gene's
  exons are those of its MANE Select transcript, else its Ensembl canonical
  transcript, else its longest coding one, so a GENCODE GTF is recommended.

### Running it

One sample, over every protein-coding gene, with 8 processes:

```
python find_domain_dups.py --bam tumor.bam --genome GRCh38.fa --gtf gencode.gtf -t 8 -o tumor
```

A tumor and its matched normal, over a list of genes (a gene name or id per line,
or a BED with the names in column 4). Every event found in the tumor is also counted
in the normal:

```
python find_domain_dups.py --bam tumor.bam --normal normal.bam --genome GRCh38.fa --gtf gencode.gtf \
    --genes cancer_genes.txt --sample_name T1 --normal_name N1 -t 8 -o T1
```

The events of an earlier run, counted in another sample, such as each library of a
panel of normals, to tell tumor-specific events from recurrent artifacts:

```
python find_domain_dups.py --bam normal2.bam --events T1.events.tsv --genome GRCh38.fa --gtf gencode.gtf -o normal2.T1
```

### Outputs

- `PREFIX.events.tsv`: one row per event, with:
  - the gene, transcript, event (`duplication` or `inversion`), exons, and copies (of
    a duplication, at least 2; 1 for an inversion)
  - the genomic span and length
  - whether it is coding and in frame, the duplicated residues, and the protein length
  - per sample (`sample_`, and `normal_` with `--normal`) the supporting reads, the
    unflanked reads (exons all inside the event, as from a circular RNA), the
    reference reads and the fraction
- `PREFIX.reads.tsv`: each read with an event and its exon order (exon numbers of the
  transcript, `r` for reverse complement)
- `PREFIX.proteins.fa`: each gene's reference protein and the protein of each coding
  duplication

### Options

| Option | Default | Meaning |
|---|---|---|
| `--min_reads` | 3 | supporting reads in `--bam` to report an event |
| `--min_insertion` | 30 | an insertion at least this long marks a read to search |
| `--min_clip` | 30 | as does a soft clip at least this long |
| `--min_exon` | 12 | shorter exons are not searched for |
| `--max_error` | 0.05 | edits per base allowed in an exon match |
| `--max_copies` | 6 | copies of an exon searched for in a read |
| `-t`, `--threads` | 1 | genes processed at once |

Run `python find_domain_dups.py --help` for the full description.

## Models of a gene's events, and the reads realigned to them (`build_dup_models.py`)

In the published analysis, the reads of a duplicated gene were realigned to transcript
models with the duplication, and variants (such as resistance mutations) were called
and phased on the duplicated copies from those alignments. `find_domain_dups.py` finds
events without this step. `build_dup_models.py` does it for one gene of its
`events.tsv`, with a model of every event found in the gene:

- **Models**: the gene's normal transcript (exons 1..N, the transcript
  `find_domain_dups.py` used) and one model per event of the gene in the
  `events.tsv`: a duplication of exons i..j in c copies is 1..j, then i..j c - 1 more
  times, then j+1..N (so a gene with a duplication and a triplication of the same
  exons gets both models), and an inversion is 1..i-1, j..i reverse complemented,
  j+1..N. The models are exactly the events the search reported, so an event to be
  modelled has to be found by `find_domain_dups.py` first.
- **Reads**: every read with an alignment over any of the events' exons, in each
  `--bam`, not only those `find_domain_dups.py` found the events in, including reads
  whose primary alignment is elsewhere in the gene. Each read's sequence and base
  qualities are taken from its primary alignment, with its haplotag (HP, PS) tags,
  which minimap2 copies to the new alignment.
- **Realignment**: minimap2 (`--secondary=no --MD`) to the models, sorted and indexed.
- **Counts**: each event is counted by the junctions reads cross, whatever the model
  minimap2 chose and its mapping quality, as a read crossing one copy-to-copy junction
  fits the duplication and triplication models alike. A duplication of c copies counts
  the reads crossing its exon j to exon i junction at least c - 1 times, so showing at
  least c copies, as `find_domain_dups.py` counts them. An inversion counts the reads
  crossing into or out of its reverse complemented exons. Normal reads are aligned to
  the normal model, with a mapping quality of at least `--min_mapq`, and cross the
  junction the event replaces. The fraction is event / (event + normal) reads.

Requirements: as `find_domain_dups.py` (it imports it from the same directory), with
minimap2 and samtools in PATH.

A tumor and its normal, every event `find_domain_dups.py` found in EGFR:

```
python build_dup_models.py --events T1.events.tsv --gene EGFR \
    --bam T1.bam --bam N1.bam --sample_names T1,N1 \
    --genome GRCh38.fa --gtf gencode.gtf -t 8 -o T1.EGFR
```

Realignment finds more event reads than the search, as it counts every
read crossing the junctions, not only those with a large insertion, clip or split
alignment.

Outputs:

- `PREFIX.models.fa` (and `.fai`): the models, named `GENE_normal`,
  `GENE_duplication_18-25`, `GENE_duplication_18-25_x3`, `GENE_inversion_5-6` and so on
- `PREFIX.models.exons.bed`: each model's exons in model coordinates, named by exon
  number, with `.2`, `.3` for further copies and `r` for reverse complement. Load it in
  IGV with the models FASTA and the BAMs
- `PREFIX.<sample>.bam` (and `.bai`): each sample's reads aligned to the models, with
  their haplotags, for calling variants with `PREFIX.models.fa` as the reference and
  phasing them
- `PREFIX.events.tsv`: per sample and event, its model, the event and normal reads,
  and the fraction
- `PREFIX.models.tsv`: per sample and model, the reads aligned to it, and those of at
  least `--min_mapq`
- `PREFIX.reads.tsv`: per read, its sample, model, mapping quality, the junctions it
  crosses (as exon pairs, `25-18` for a step back) and the events it supports

Run `python build_dup_models.py --help` for every option.
