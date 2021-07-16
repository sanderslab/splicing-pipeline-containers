from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

@pipeline_task
def run_exon_expression(task_context):
    '''
    Runs HTSeq/DEXseq for getting exon read counts
    '''
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix = my_line[0]
    bam_uri = my_line[1]
    gene_name = my_line[2]
    region = my_line[3]

    out_path = task_context.path_for_volume_file(f'{prefix}_{gene_name}.txt')

    bam_path = task_context.download_s3_object(bam_uri)
    sub_bam_path = task_context.path_for_volume_file(f'{prefix}_{gene_name}.bam')
    gff_path = task_context.get_asset(f'gencode.v31.annotation.{gene_name}.gff')

    sys.stdout.flush()
    print(f':::INDEXING {prefix}:::')
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
    print(f':::SUBSETTINB {prefix} TO {region}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'samtools view -h {bam_path} {region} > {sub_bam_path}'
        ),
        mappings={
            'bam_path' : bam_path,
            'region' : region,
            'sub_bam_path' : sub_bam_path
        }
    )

    sys.stdout.flush()
    print(f':::GETTING COUNTS FOR {prefix}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'python3.6 /home/dexseq_count.py '
            '{gff_path} '
            '{sub_bam_path} '
            '{out_path}'
        ),
        mappings={
            'gff_path' : gff_path,
            'sub_bam_path' : sub_bam_path,
            'out_path' : out_path
        }
    )

    sys.stdout.flush()
    print(f':::FINISHED; UPLOADING {out_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('exon_expression', 'txt', out_path)
    task_context.put_result_for_stage('exon_expression', 'bam', sub_bam_path)
    sys.stdout.flush()

if __name__ == '__main__':
    run_exon_expression()
