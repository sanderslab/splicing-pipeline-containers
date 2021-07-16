from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

def tester(task_context):
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    uri = my_line[0]
    path = task_context.download_s3_object(uri)
    task_context.put_result_for_stage('gene_expression', 'txt', path)


@pipeline_task
def run_gene_expression(task_context):
    '''
    Runs HTSeq/DEXseq for getting exon read counts
    '''
    # tester(task_context)
    # return
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix, bam_uri, ended_type = my_line

    if ended_type == 'STR':
        strandedness = 'yes'
    elif ended_type == 'UNSTR':
        strandedness = 'no'
    elif ended_type == 'REV':
        strandedness = 'reverse'
    else:
        raise ValueError(f'Unrecognized end type {ended_type}')

    out_path = task_context.path_for_volume_file(f'{prefix}_gene_exp.txt')

    sorted_bam_path = task_context.download_s3_object(bam_uri)
    # bam_path = task_context.get_result_for_stage('star', 'bam').get_array_result().download_to_volume()
    # sorted_bam_path = bam_path.replace('.bam', '.sorted.bam')
    gtf_path = task_context.get_asset(f'gencode.v31.annotation.gtf')
    # gtf_path = task_context.get_asset(f'gencode.vM25.chr_patch_hapl_scaff.annotation.gtf')
    # gtf_path = task_context.get_asset(f'rheMac8.ensGene.gtf')
    sys.stdout.flush()
    print(f':::INDEXING {prefix}:::')
    sys.stdout.flush()

    sys.stdout.flush()
    print(f':::GETTING COUNTS FOR {prefix}:::')
    sys.stdout.flush()

    task_context.run_template(
        #time htseq-count bam gtf -r pos -s yes > P11_HET_F_S27Aligned.htseq.out
        template=(
            'htseq-count '
            '{bam_path} '
            '{gtf_path} '
            '-s no '
            '-r pos '
            '> {out_path}'
        ),
        mappings={
            'strandedness' : strandedness,
            'gtf_path' : gtf_path,
            'bam_path' : sorted_bam_path,
            'out_path' : out_path
        }
    )

    sys.stdout.flush()
    print(f':::FINISHED; UPLOADING {out_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('gene_expression', 'txt', out_path)
    # task_context.put_result_for_stage('sorted_bam', 'bam', sorted_bam_path)
    # task_context.put_result_for_stage('sorted_bam', 'bai', f'{sorted_bam_path}.bai')

    sys.stdout.flush()

if __name__ == '__main__':
    run_gene_expression()
