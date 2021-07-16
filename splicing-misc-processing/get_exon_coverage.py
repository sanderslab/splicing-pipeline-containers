'''
Using bedtools, get coverage of each individual exon
Right now this script is very specific to brainvar and SCN2A
'''
from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

LOCUS = 'chr2:165210285-165416691'

@pipeline_task
def extract(task_context):
    '''
    Manifest for this stage will be the same as olego manifest
    '''
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix = my_line[0]
    bam_uri = my_line[1]

    bam_path = task_context.download_s3_object(bam_uri)
    bed_path = task_context.get_asset('SCN2A_exons.bed')
    out_path = task_context.path_for_volume_file(f'{prefix}_{LOCUS}_counts.bed')

    sys.stdout.flush()
    print(f':::INDEXING {bam_path}:::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'samtools index {bam_path}'
        ),
        mappings={
            'bam_path' : bam_path
        }
    )

    sys.stdout.flush()
    print(f':::EXTRACTING {LOCUS}:::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'samtools view -b {bam_path} {locus} -o {locus_bam}'
        ),
        mappings={
            'bam_path' : bam_path,
            'locus' : LOCUS,
            'locus_bam' : f'{prefix}_{LOCUS}.bam'
        }
    )

    sys.stdout.flush()
    print(f':::RUNNING BEDTOOLS:::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            ' bedtools intersect -b {locus_bam} -a {bed} -bed -c > {out_path}'
        ),
        mappings={
            'locus_bam' : f'{prefix}_{LOCUS}.bam',
            'bed' : bed_path,
            'out_path' : out_path
        }
    )

    sys.stdout.flush()
    print(':::UPLOADING RESULTS:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('exon_cov', 'bed_result', out_path)
    task_context.put_result_for_stage('exon_cov', 'bam', f'{prefix}_{LOCUS}.bam')
    task_context.put_result_for_stage('exon_cov', 'bam', f'{prefix}_{LOCUS}.bam.bai')
    sys.stdout.flush()

if __name__ == '__main__':
    extract()


