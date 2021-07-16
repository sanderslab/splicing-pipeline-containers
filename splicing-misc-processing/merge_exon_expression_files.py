from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import csv
import itertools
import os
import sys

def __get_exon_positions(gff_path):
    '''
    Returns a dict of exon positions
    '''
    exon_dict = dict()
    with open(gff_path, 'r') as gff_in:
        for line in gff_in:
            chrom, _, locus_type, start, end, _, strand, _, info = [i.strip() for i in line.split('\t')]
            if locus_type == 'exonic_part':
                locus = f'{chrom}:{start}-{end}'
                gene_id = [i.strip() for i in info.split(';')][0].split('"')[1]
                exon_num = [i.strip() for i in info.split(';')][-1].split('"')[1]
                exon_id = f'{gene_id}:{exon_num}'
                exon_dict[exon_id] = locus

    return exon_dict

def __get_gene_positions(gtf_path):
    '''
    For HTSeq runs
    '''
    gene_dict = dict()
    with open(gtf_path, 'r') as gtf_in:
        for line in gtf_in:
            if not line.startswith('#'):
                chrom, _, locus_type, start, end, _, strand, _, info = [i.strip() for i in line.split('\t')]
                if locus_type == 'gene':
                    locus = f'{chrom}:{start}-{end}'
                    gene_id = [i.strip() for i in info.split(';')][0].split('"')[1]
                    gene_dict[gene_id] = locus

    return gene_dict

def __make_merge_file(in_files_paths, out_file):
    '''
    Make a tsv of exon read counts per each sample
    '''

    samp_exon_exp = dict()
    with open(out_file, 'w') as out_f:
        w = csv.writer(out_f, delimiter='\t')
        for file_path in in_files_paths:
            samp = file_path.split('/')[-1].split('.')[0]
            with open(file_path, 'r') as in_f:
                for line in in_f:
                    if not line.startswith('_'):
                        exon, count = [i.strip() for i in line.split('\t')]
                        if samp not in samp_exon_exp:
                            samp_exon_exp[samp] = {exon: count}
                        else:
                            samp_exon_exp[samp][exon] = count

        exons = sorted(list(samp_exon_exp.values())[0].keys())
        w.writerow(['#Sample'] + exons)
        for samp in samp_exon_exp:
            w.writerow([samp] + [samp_exon_exp[samp][exon] for exon in exons])

def __annotate_exons(in_file, exon_dict, out_file):
    '''
    Annotate header line of file made in make_merge_file with 
    exon dict from get_exon positions
    '''
    locus_header = list()
    with open(out_file, 'w') as anno_f:
        w = csv.writer(anno_f, delimiter='\t')
        with open(in_file, 'r') as in_f:
            for line in in_f:
                if line.startswith('#'):
                    header = [i.strip() for i in line.split('\t')]
                    for item in header:
                        if item.startswith('ENS'):
                            locus = exon_dict[item]
                            locus_header.append(locus)
                    w.writerow(['#Sample'] + locus_header)
                # w.writerow(line.split('\t'))
                anno_f.write(line)



@pipeline_task
def merge_exon_expression_files(task_context):
    '''
    Runs HTSeq/DEXseq for getting exon read counts
    '''
    # my_line = tuple(task_context.get_my_input_line().strip().split('\t'))

    txt_objs = task_context.get_result_for_stage('exon_expression', 'txt')
    txt_paths = [txt.download_to_volume() for txt in txt_objs]

    # gff_path = task_context.get_asset(f'gencode.v31.annotation.gff')
    gtf_path = task_context.get_asset(f'gencode.v31.annotation.gtf')
    # gff_path = task_context.get_asset(f'gencode.v31.annotation_SYNGAP_exon31_alt.gff')
    # gff_path = task_context.get_asset(f'gencode.vM25.chr_patch_hapl_scaff.annotation.gff')
    # gtf_path = task_context.get_asset(f'gencode.vM25.chr_patch_hapl_scaff.annotation.gtf')
    # gff_path = task_context.get_asset(f'rheMac8.ensGene.gff')

    merged_path = task_context.path_for_volume_file('cohort_gene_exp.txt')
    anno_path = task_context.path_for_volume_file('cohort_gene_exp.anno.txt')

    sys.stdout.flush()
    print(f':::MERGING:::')
    sys.stdout.flush()

    # exon_dict = __get_exon_positions(gff_path)
    gene_dict = __get_gene_positions(gtf_path)
    __make_merge_file(txt_paths, merged_path)
    # __annotate_exons(merged_path, exon_dict, anno_path)
    __annotate_exons(merged_path, gene_dict, anno_path)

    sys.stdout.flush()
    print(f':::FINISHED; UPLOADING {anno_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('merged_gene_expression', 'txt', merged_path)
    task_context.put_result_for_stage('merged_gene_expression', 'anno', anno_path)
    sys.stdout.flush()

if __name__ == '__main__':
    merge_exon_expression_files()
