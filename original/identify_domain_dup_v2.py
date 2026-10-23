import os
# import subprocess
import pipettor
import pysam
import edlib

BEDTOOLS_PATH = '/private/groups/brookslab/cafelton/bin/bedtools2/bin/bedtools'
min_read_support = 3

def get_reverse_complement(seq):
    compbase = {'A': 'T', 'T': 'A', 'C': 'G', 'G': 'C', 'N': 'N',
                'R': 'Y', 'Y': 'R', 'K': 'M', 'M': 'K', 'S': 'S',
                'W': 'W', 'B': 'V', 'V': 'B', 'D': 'H', 'H': 'D'}
    seq = seq.upper()
    new_seq = []
    for base in seq:
        new_seq.append(compbase[base])
    return ''.join(new_seq[::-1])

ensgtoname = {}
for line in open('/private/groups/brookslab/cafelton/bin/gencode38.ensg.hugo.tsv'):
    line = line.rstrip().split('\t')
    ensgtoname[line[0].split('.')[0]] = line[1]

gene_conv = {'MYCL1':'MYCL', 'GLTSCR1':'BICRA', 'GLTSCR2':'NOP53', 'MRE11A':'MRE11', 'FAM175A':'ABRAXAS1', 'HIST1H3B':'H3C2', 'H3F3B':'H3-3B', 'H3F3A':'H3-3A', 'C11ORF95':'ZFTA'}
onco_genes = set()
for line in open('oncoplex_v6_gene_list.csv'):
    line = line.rstrip().lstrip('\ufeff')
    if line in gene_conv:
        line = gene_conv[line]
    onco_genes.add(line)

for line in open('cancerGeneList.tsv'):
    line = line.split('\t')
    onco_genes.add(line[0])


batch_to_name = {'Berger_2509_PB1283': 'D_2509', 'Berger_2507_PB1264': 'D_2507', 'Berger_patientdata_110124': 'UW_2411', 'Berger_042225': 'UW_2504', 'Berger_2508_PB1281': 'D_2508', 'Castle_PacBio_0925': 'D_C2509', 'Berger_cellline_1124': 'D_2411'}

all_samples = set()
for line in open('/private/groups/brookslab/cafelton/Berger_collab/1025_data_processing/Berger_raw_data_manifest_101325_updatedIDs_withCLfix.tsv'):
    if line[:4] != 'Data':
        line = line.rstrip().split('\t')
        folder = line[1]
        file = line[2]
        newfile = line[3]
        batch = batch_to_name[folder.split('/')[-2]]
        sample = '_'.join(file.split('.')[0].split('_')[:-1])
        newsample = '_'.join(newfile.split('.')[0].split('_')[:-1])
        sample = batch + '_' + sample
        newsample = batch + '_' + newsample

        newsample = tuple(newsample.split('_'))
        all_samples.add(newsample)

print('loading annotation')
gtf_file = '/private/groups/brookslab/reference_annotations/gencode.v49.basic.annotation.gtf'
gene_to_exons = {}
gene_to_bounds = {}
for line in open(gtf_file):
    if line[0] != '#':
        line = line.split('\t')
        if line[2] == 'gene':
            gene = line[-1].split('gene_id "')[1].split('"')[0]
            gene_to_bounds[gene] = (line[0], line[6], int(line[3])-1, int(line[4]))
        if line[2] == 'exon' and 'MANE_Select' in line[-1]:
            gene = line[-1].split('gene_id "')[1].split('"')[0]
            if gene not in gene_to_exons:
                gene_to_exons[gene] = [line[0], line[6], []]
            gene_to_exons[gene][-1].append((int(line[3])-1, int(line[4])))

gene_bounds_file = gtf_file.split('/')[-1].split('.gtf')[0] + '.gene_bounds.bed'

with open(gene_bounds_file, 'w') as out:
    for gene in gene_to_bounds:
        bounds = gene_to_bounds[gene]
        out.write('\t'.join([str(x) for x in [bounds[0], bounds[2], bounds[3], gene, '.', bounds[1]]]) + '\n')

all_samples = sorted(list(all_samples))

#identify what gene the insertion belongs to
# get the exon sequence for that gene
# extract the reads in the gene locus from the bam

genome = pysam.FastaFile('/private/groups/brookslab/reference_sequence/GRCh38.p14.genome.fa')

out_final = open(f'funcotator_0226/all.dup.final.050626.txt', 'a')

c = 0
for s in all_samples:
    seqcenter, batch, seqtype, source_type, tissue_type, tumor_state, sample_name = s
    sample_id = '_'.join(s)
    if seqtype == 'PB':
        temp_prefix = f'{sample_id}.temp'
        out_prefix = f'funcotator_0226/{sample_id}'
        
        if ((tumor_state == 'T' and source_type == 'Ti') or source_type == 'CL' or tumor_state == 'TR'): #and sample_name in {'PC9', 'MORE19post', 'ULC1627T'}: #and sample_name == 'ULC0563T': #'PC9':#'ULC0513T': #and sample_name == 'MORE19post': #and sample_name == 'ULC0513N': #
            c += 1
            if c < 24:
                continue
            large_insert_file = f'{out_prefix}.indels.exonOverlap.largeInsertions.bed'
            # print(large_insert_file)
            bam_file = f'genome_aligned/{sample_id}.genomealigned.withMD.bam'
            genomic_chim_file = f'/scratch/cafelton/berger/genome_chim/{sample_id}.genomealigned.chim.bam'
            # transcriptomic_chim_file = f'/scratch/cafelton/berger/tome_chim/{sample_id}.gencode49basicaligned.chim.sorted.bam'

            # pipettor.run([(f'{BEDTOOLS_PATH}', 'bamtobed', '-i', genomic_chim_file)], stdout=f'{temp_prefix}.genomicchim.bed')
            # pipettor.run([(f'{BEDTOOLS_PATH}', 'intersect', '-a', f'{temp_prefix}.genomicchim.bed', '-b', gene_bounds_file, '-wo')], stdout=f'{temp_prefix}.genomicchim.geneoverlap.txt')
            align_to_genes = {}
            for line in open(f'{temp_prefix}.genomicchim.geneoverlap.txt'):
                line = line.split('\t')
                gene = line[9]
                read = tuple(line[:4])
                if read not in align_to_genes:
                    align_to_genes[read] = []
                align_to_genes[read].append(gene)
            read_to_genes = {}
            for chrom, start, end, read in align_to_genes:
                if read not in read_to_genes:
                    read_to_genes[read] = []
                read_to_genes[read].append(align_to_genes[(chrom, start, end, read)])
            gene_to_reads = {}
            for read in read_to_genes:
                pot_combos = [[]]
                for g_list in read_to_genes[read]:
                    new_pc = []
                    for pc in pot_combos:
                        for g in g_list:
                            new_pc.append(pc + [g])
                    pot_combos = new_pc
                for comb in pot_combos:
                    if len(comb) > 1 and len(set(comb)) == 1:
                        gene = comb[0]
                        if gene not in gene_to_reads:
                            gene_to_reads[gene] = 0
                        gene_to_reads[gene] += 1
            prelim_genes = {}
            for gene in gene_to_reads:
                if gene_to_reads[gene] >= 10:
                    prelim_genes[gene] = gene_to_reads[gene]
            # print('prelim', len(prelim_genes), list(prelim_genes)[:10])
                
            if os.path.getsize(large_insert_file) > 0:
                print(c, sample_name, 'identifying genes with insertions')
                # pipettor.run([(f'{BEDTOOLS_PATH}', 'intersect', '-a', large_insert_file, '-b', gene_bounds_file, '-wo')], stdout=f'{temp_prefix}.largeInserts.withGenes.txt')
                for line in open(f'{temp_prefix}.largeInserts.withGenes.txt'):
                    line = line.split('\t')
                    gene = line[7]
                    rcount = int(line[3].split('_')[-1].split('/')[0])
                    if gene not in prelim_genes:
                        prelim_genes[gene] = rcount
                    else:
                        prelim_genes[gene] += rcount
                
            insert_genes = {}
            for gene in prelim_genes:
                if gene not in insert_genes and gene in gene_to_exons:
                    chrom, gstrand, gstart, gend = gene_to_bounds[gene]
                    exon_seq = []
                    for estart, eend in sorted(gene_to_exons[gene][-1]):
                        exon_seq.append(genome.fetch(chrom, estart, eend))
                    if gstrand == '-':
                        exon_seq = exon_seq[::-1]
                        exon_seq = [get_reverse_complement(x) for x in exon_seq]
                    full_exon_seq = {}
                    for exon_index, ref_exon in enumerate(exon_seq):
                        full_exon_seq['fwd_' + str(exon_index + 1)] = ref_exon
                        full_exon_seq['rev_' + str(exon_index + 1)] = get_reverse_complement(ref_exon)
                    
                    # print(chrom, gstrand, gstart, gend, [len(x) for x in exon_seq])
                    insert_genes[gene] = (chrom, gstrand, gstart, gend, full_exon_seq)
                
            my_genes = []
            for gene in insert_genes:
                chrom, gstrand, gstart, gend, exon_seq = insert_genes[gene]
                
                tot_reads = 0
                for a in pysam.AlignmentFile(bam_file, 'rb').fetch(chrom, gstart, gend):
                    if not a.is_secondary:
                        tot_reads += 1
                gname = ensgtoname[gene.split('.')[0]] if gene.split('.')[0] in ensgtoname else gene.split('.')[0]
                
                if prelim_genes[gene]/tot_reads > 0.01 and gname in onco_genes:
                    # print(gname, prelim_genes[gene], tot_reads, prelim_genes[gene]/tot_reads)
                    my_genes.append(gene)
                    # print(gene, [len(exon_seq[x]) for x in exon_seq])
            # print('sufficient reads + onco:', len(my_genes))
            # break
            for gene in my_genes:
                chrom, gstrand, gstart, gend, exon_seq = insert_genes[gene]
                gname = ensgtoname[gene.split('.')[0]] if gene.split('.')[0] in ensgtoname else gene.split('.')[0]
                print(c, sample_name, 'processing reads for', gname)
                if gname in ['AFDN', 'ATRX', 'PSIP1', 'CCAR1', 'TOP1']:
                    continue
                # if gname not in ['AFDN', 'ATRX', 'TOP1', 'MAPKAP1']:
                #     continue


                all_exon_orders = {}
                for a in pysam.AlignmentFile(bam_file, 'rb').fetch(chrom, gstart, gend):
                    if not a.is_secondary:
                        rstrand = '-' if a.is_reverse else '+'
                        if rstrand == gstrand:
                            rseq = a.get_forward_sequence()
                            og_rseq = rseq
                            all_exon_locs = []
                            for exon_index, ref_exon in exon_seq.items():
                                if len(ref_exon) >= 20: ##if exon is too short, edlib can actually kill the computer trying to match the short seq
                                    my_aligns = None
                                    while (my_aligns == None or (my_aligns['editDistance'] != -1 and my_aligns['editDistance'] < 5)): #and d < 5:
                                        my_aligns = edlib.align(ref_exon, rseq, mode='HW', task='locations', k=5)
                                        # print(gname, my_aligns)
                                        rseq = list(rseq)
                                        # if len(my_aligns['locations']) > 1:
                                        #     print(a.query_name, exon_index, my_aligns['locations'])
                                        added = []
                                        for i, (s, e) in enumerate(my_aligns['locations']):
                                            if s == None:
                                                s = 0
                                            if i == 0 or not min(e, my_aligns['locations'][0][1]) - max(s, my_aligns['locations'][0][0]) > 0:
                                                all_exon_locs.append((s, e, exon_index)) ##adding 1 for 1-based number
                                                rseq[s:e+1] = ['N']* ((e+1)-s)
                                                added.append((s, e))
                                        # if len(added) > 1:
                                        #     print(a.query_name, exon_index, my_aligns['locations'], added)
                                        rseq = ''.join(rseq)
                            all_exon_locs.sort()
                            final_locs = []
                            lastend = 0
                            for s, e, i in all_exon_locs:
                                if s - lastend > 5:
                                    final_locs.append((-1, -1,'?'))
                                final_locs.append((s, e, i))
                                lastend=e
                            for i in range(1, len(final_locs)-2):
                                if final_locs[i][2] == '?':
                                    ##if missing exon is clearly sequentially between adjacent exons, correct to exon number
                                    prev_exon, next_exon = final_locs[i-1][2], final_locs[i+1][2]
                                    if prev_exon.split('_')[0] == next_exon.split('_')[0]:
                                        if prev_exon.split('_')[0] == 'fwd' and int(prev_exon.split('_')[1]) + 2 == int(next_exon.split('_')[1]):
                                            final_locs[i] = (final_locs[i][0], final_locs[i][1], 'fwd_' + str(int(prev_exon.split('_')[1]) + 1))
                                        if prev_exon.split('_')[0] == 'rev' and int(prev_exon.split('_')[1]) - 2 == int(next_exon.split('_')[1]):
                                            final_locs[i] = (final_locs[i][0], final_locs[i][1], 'rev_' + str(int(prev_exon.split('_')[1]) - 1))
                            final_locs_short = [x[2] for x in final_locs]
                            if len(final_locs_short) >= 2:
                                # if final_locs_short[0] == '?':
                                #     final_locs_short = final_locs_short[1:]
                                # if final_locs_short[-1] == '?':
                                #     final_locs_short = final_locs_short[:-1]
                                final_locs_short = [x for x in final_locs_short if x != '?']
                                if len(final_locs_short) >= 2:
                                    final_locs_short = tuple(final_locs_short)
                                    if final_locs_short not in all_exon_orders:
                                        all_exon_orders[final_locs_short] = 0
                                    all_exon_orders[final_locs_short] += 1
                
                # for eorder in sorted(list(all_exon_orders.keys()), key=lambda x:len(x), reverse=True):
                #     print(eorder, all_exon_orders[eorder])
                #     # if all_exon_orders[eorder] >= min_read_support:
                #     #     print(eorder, all_exon_orders[eorder])
                # break

                confident_long_eorders = {}
                for eorder in sorted(list(all_exon_orders.keys()), key=lambda x:len(x), reverse=True):
                    if all_exon_orders[eorder] >= min_read_support:
                        seen = False
                        for e2 in confident_long_eorders:
                            if str(eorder)[1:-1] in str(e2):
                                seen = True
                        if not seen:
                            confident_long_eorders[eorder] = all_exon_orders[eorder]

                dup_exons = set()
                has_dup = False
                for e in confident_long_eorders:
                    for x in e:
                        if e.count(x) > 1 or x[:3] == 'rev':
                            has_dup = True
                            dup_exons.add(x.split('_')[-1])
                # print(dup_exons)
                final_eorders = [tuple([f'fwd_{x}' for x in range(1, len(gene_to_exons[gene][-1]) + 1)])] #initialize with normal MANE transcript
                for e in confident_long_eorders:
                    if len(set([x.split('_')[-1] for x in e]) & dup_exons) > 0: # only keep models that cover duplicated exons
                        new_order = []
                        # add exons to start and end to make full transcript
                        first_exon, last_exon = int(e[0].split('_')[-1]), int(e[-1].split('_')[-1])
                        if first_exon > 1:
                            new_order.extend([f'fwd_{x}' for x in range(1, first_exon)])
                        new_order.extend(e)
                        if last_exon < len(gene_to_exons[gene][-1]):
                            new_order.extend([f'fwd_{x}' for x in range(last_exon+1, len(gene_to_exons[gene][-1]) + 1)])
                        
                        #add exons to remove exon skipping
                        new_order2 = [new_order[0]]
                        for exon in new_order[1:]:
                            last_pos = int(new_order2[-1].split('_')[-1])
                            this_pos = int(exon.split('_')[-1])
                            if this_pos - last_pos > 1:
                                for i in range(last_pos+1, this_pos):
                                    new_order2.append(f'fwd_{i}')
                            elif this_pos - last_pos <= 0:
                                exons_between = []
                                for other_eorder in confident_long_eorders:
                                    if f'fwd_{last_pos}' in other_eorder and f'fwd_{this_pos}' in other_eorder:
                                        first_index = -1
                                        for i in range(len(other_eorder)):
                                            if other_eorder[i] == f'fwd_{this_pos}' and first_index != -1:
                                                these_exons = other_eorder[first_index + 1:i]
                                                if len(these_exons) > len(exons_between):
                                                    exons_between = these_exons
                                                first_index = -1
                                            if other_eorder[i] == f'fwd_{last_pos}':
                                                first_index = i
                                            
                                    # try:
                                    #     first_index = other_eorder.index(f'fwd_{last_pos}')
                                    # except ValueError:
                                    #     first_index = -1
                                    # try:
                                    #     last_index = other_eorder.index(f'fwd_{this_pos}')
                                    # except ValueError:
                                    #     last_index = -1
                                    # if first_index != -1 and last_index != -1 and first_index < last_index:
                                    #     these_exons = other_eorder[first_index+1:last_index]
                                    #     if len(these_exons) > len(exons_between):
                                    #         exons_between = these_exons
                                new_order2.extend(exons_between)



                                # pot_add_exon = last_pos + 1
                                # while str(pot_add_exon) in dup_exons:
                                #     new_order2.append(f'fwd_{pot_add_exon}')
                                #     pot_add_exon += 1
                                # pot_add_exon = this_pos - 1
                                # exons_to_add = []
                                # while str(pot_add_exon) in dup_exons:
                                #     exons_to_add.insert(0, f'fwd_{pot_add_exon}')
                                #     pot_add_exon -= 1
                                # new_order2.extend(exons_to_add)
                            new_order2.append(exon)
                        
                        
                        new_order = tuple(new_order2)
                        if new_order not in final_eorders:
                            final_eorders.append(tuple(new_order))
                            print(gene, ', '.join(new_order))


                if has_dup:

                    dup_ref = f'{out_prefix}.{gname}.dup.transcriptome-2.fa'
                    dup_annot = f'{out_prefix}.{gname}.dup.transcriptome.annotexons-2.bed'
                    dup_aligned = f'{out_prefix}.{gname}.dup.transcriptome.aligned-2.bam'
                    dup_normal_aligned = None
                    bam_file_normal = None
                    if tumor_state == 'T' and source_type == 'Ti':
                        normal_id = (seqcenter, batch, seqtype, source_type, tissue_type, 'N', sample_name[:-1] + 'N')
                        normal_id = '_'.join(normal_id)
                        bam_file_normal = f'genome_aligned/{normal_id}.genomealigned.withMD.bam'
                        dup_normal_aligned = f'funcotator_0226/{sample_id}.{gname}.dup.transcriptome.normal.aligned-2.bam'
                    print(c, sample_name, 'generating and aligning to domain dup models for', gname)
                    
                    model_to_exon_pos = {}
                    with open(dup_ref, 'w') as faout, open(dup_annot, 'w') as annotout:
                        # isoid = 1
                        for e in final_eorders:
                            dup_groups = []
                            has_dup = False
                            for x in e:
                                if e.count(x) > 1 or x[:3] == 'rev':
                                    has_dup = True
                            if has_dup:
                                dup_name = ['dup']
                            else:
                                dup_name = ['none']
                            this_group = [e[0], e[0]]
                            for x in e[1:]:
                                if x.split('_')[0] == this_group[-1].split('_')[0]:
                                    if x.split('_')[0] == 'fwd' and int(x.split('_')[-1]) == int(this_group[-1].split('_')[-1]) + 1:
                                        this_group[-1] = x
                                    elif x.split('_')[0] == 'rev' and int(x.split('_')[-1]) == int(this_group[-1].split('_')[-1]) - 1:
                                        this_group[-1] = x
                                    else:
                                        dup_groups.append(tuple(this_group))
                                        this_group = [x, x]
                                else:
                                    dup_groups.append(tuple(this_group))
                                    this_group = [x, x]
                            dup_groups.append(tuple(this_group))
                        
                            for g in dup_groups:
                                if g[0] == g[1]:
                                    dup_name.append(f"{':'.join(g[0].split('_'))}")
                                else:
                                    dup_name.append(f"{g[0].split('_')[0]}:{'-'.join([x.split('_')[-1] for x in g])}")
                            dup_name = '_'.join(dup_name)
                            # print(e)
                            # print(dup_name)

                            faout.write('>' + gname + '_' + dup_name + '\n')
                            faout.write(''.join([exon_seq[x] for x in e if x != '?']) + '\n')
                            currtot = 0
                            model_to_exon_pos[gname + '_' + dup_name] = []
                            for x in e:
                                if x != '?':
                                    model_to_exon_pos[gname + '_' + dup_name].append((currtot, currtot + len(exon_seq[x]), x))
                                    outline = [gname + '_' + dup_name, currtot, currtot + len(exon_seq[x]), x]
                                    annotout.write('\t'.join([str(x) for x in outline]) + '\n')
                                    currtot += len(exon_seq[x])
                            # isoid += 1
                    pipettor.run([('samtools', 'faidx', dup_ref)])
                    get_region_cmd = ('samtools', 'view', '-hb', bam_file, f'{chrom}:{gstart}-{gend}')
                    get_seq_cmd = ('samtools', 'fasta', '-')
                    minimap_cmd = ('minimap2', '-a', '--secondary=no', dup_ref, '-')
                    conv_bam_cmd = ('samtools', 'view', '-hb', '-')
                    sort_cmd = ('samtools', 'sort', '-o', dup_aligned, '-')
                    pipettor.run([get_region_cmd, get_seq_cmd, minimap_cmd, conv_bam_cmd, sort_cmd])
                    pysam.index(dup_aligned)

                    if dup_normal_aligned is not None:
                        get_region_cmd = ('samtools', 'view', '-hb', bam_file_normal, f'{chrom}:{gstart}-{gend}')
                        sort_cmd = ('samtools', 'sort', '-o', dup_normal_aligned, '-')
                        pipettor.run([get_region_cmd, get_seq_cmd, minimap_cmd, conv_bam_cmd, sort_cmd])
                        pysam.index(dup_normal_aligned)

                    dup_to_reads = {}
                    infile = pysam.AlignmentFile(dup_aligned, 'rb')
                    outfile = pysam.AlignmentFile(dup_aligned.split('.bam')[0] + '.filtered.bam', 'wb', template=infile)
                    exon_cov_to_reads = {}
                    for a in infile:
                        if a.is_mapped and not a.is_secondary and not a.is_supplementary: #and a.mapping_quality > 1:
                            tot_insert_len = a.get_cigar_stats()[0][1]
                            if tot_insert_len < 50:
                                dup = a.reference_name
                                covered_exons = []
                                for s, e, n in model_to_exon_pos[dup]:
                                    if a.reference_start < s and e < a.reference_end:
                                        covered_exons.append(n)
                                    elif (a.reference_start < e - ((e-s)/2) and e < a.reference_end) \
                                      or (a.reference_start < s and s + ((e-s)/2) < a.reference_end):
                                        covered_exons.append(n)
                                covered_exons = tuple(covered_exons)
                                if covered_exons not in exon_cov_to_reads:
                                    exon_cov_to_reads[covered_exons] = [0,0]
                                exon_cov_to_reads[covered_exons][1] += 1
                            
                            if a.mapping_quality > 1:
                                if dup not in dup_to_reads:
                                    dup_to_reads[dup] = [0, 0]
                                dup_to_reads[dup][1] += 1
                                outfile.write(a)
                    infile.close()
                    outfile.close()
                    pysam.index(dup_aligned.split('.bam')[0] + '.filtered.bam')

                    if dup_normal_aligned is not None:
                        infile = pysam.AlignmentFile(dup_normal_aligned, 'rb')
                        outfile = pysam.AlignmentFile(dup_normal_aligned.split('.bam')[0] + '.filtered.bam', 'wb', template=infile)
                        for a in infile:
                            if a.is_mapped and not a.is_secondary and not a.is_supplementary: #and a.mapping_quality > 1:
                                tot_insert_len = a.get_cigar_stats()[0][1]
                                if tot_insert_len < 50:
                                    dup = a.reference_name
                                    covered_exons = []
                                    for s, e, n in model_to_exon_pos[dup]:
                                        if a.reference_start < s and e < a.reference_end:
                                            covered_exons.append(n)
                                        elif (a.reference_start < e - ((e-s)/2) and e < a.reference_end) \
                                        or (a.reference_start < s and s + ((e-s)/2) < a.reference_end):
                                            covered_exons.append(n)
                                    covered_exons = tuple(covered_exons)
                                    if covered_exons not in exon_cov_to_reads:
                                        exon_cov_to_reads[covered_exons] = [0,0]
                                    exon_cov_to_reads[covered_exons][0] += 1
                                
                                if a.mapping_quality > 1:
                                    if dup not in dup_to_reads:
                                        dup_to_reads[dup] = [0, 0]
                                    dup_to_reads[dup][0] += 1
                                    outfile.write(a)
                        infile.close()
                        outfile.close()
                        pysam.index(dup_normal_aligned.split('.bam')[0] + '.filtered.bam')


                    model_tot, dup_tot = [0,0], [0,0]
                    for exons in exon_cov_to_reads:
                        # print(exon_cov_to_reads[exons], ', '.join(exons))
                        
                        if len([x for x in exons if x.split('_')[-1] in dup_exons]) > 0:
                            if len([x for x in exons if 'rev' in x]) > 0 or len(set(exons)) < len(exons): #has dup
                                dup_tot = [dup_tot[x] + exon_cov_to_reads[exons][x] for x in range(2)]
                        
                            model_tot = [model_tot[x] + exon_cov_to_reads[exons][x] for x in range(2)]
                        
                        # if len([x for x in exons if 'rev' in x]) > 0 or len(set(exons)) < len(exons): #has dup
                        #     dup_tot[0] += exon_cov_to_reads[exons]
                        #     if len([x for x in exons if x.split('_')[-1] in dup_exons]) > 0:
                        #         dup_tot[1] += exon_cov_to_reads[exons]
                        # model_tot[0] += exon_cov_to_reads[exons]
                        # if len([x for x in exons if x.split('_')[-1] in dup_exons]) > 0:
                        #     model_tot[1] += exon_cov_to_reads[exons]
                    print(dup_tot, model_tot)

                    if dup_tot[0] == 0 and dup_tot[1] >= min_read_support: #and dup_tot[1]/model_tot[1] >= 0.01:
                        print(gname, 'good')
                        out_final.write(f'{sample_name}\t{gname}\t{",".join(sorted(dup_exons))}\t{dup_tot[1]}\t{model_tot[1]}\t{round(dup_tot[1]/model_tot[1], 3)}\n')
                    
                    
                    # for dup in dup_to_reads:
                    #     if dup_to_reads[dup][1] >= min_read_support:
                    #         good_dup.append(dup)
                    #         model_tot = [model_tot[x] + dup_to_reads[dup][x] for x in range(2)]
                    #         if 'none' not in dup:
                    #             dup_tot = [dup_tot[x] + dup_to_reads[dup][x] for x in range(2)]
                    # print(dup_tot, model_tot)
                    # for dup in good_dup:    
                    #     print(dup, dup_to_reads[dup], round(dup_to_reads[dup][1]/model_tot[1], 3), round(dup_tot[1]/model_tot[1], 3))
                    
                    # if dup_tot[1] >= min_read_support and dup_tot[0] <= min_read_support and dup_tot[1]/model_tot[1] >= 0.1:
                    #     # print(dup_tot, model_tot)
                    #     for dup in good_dup:
                    #         # print(dup, dup_to_reads[dup], round(dup_to_reads[dup][1]/model_tot[1], 3), round(dup_tot[1]/model_tot[1], 3))
                    #         out_final.write(f'{sample_name}\t{gname}\t{dup}\t{dup_to_reads[dup][1]}\t{round(dup_to_reads[dup][1]/model_tot[1], 3)}\n')
            # out_final.close()



