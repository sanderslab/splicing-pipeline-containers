from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

def tester(task_context):
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    uri = my_line[0]
    path = task_context.download_s3_object(uri)
    task_context.put_result_for_stage('exon_expression', 'txt', path)


@pipeline_task
def run_exon_expression(task_context):
    '''
    Runs HTSeq/DEXseq for getting exon read counts
    '''
    # tester(task_context)
    # return
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    # prefix, bam_uri, ended_type = my_line
    prefix  = my_line[0]
    read_uris = my_line[1:-1]
    ended_type = my_line[-1] # SE vs PE

    if ended_type == 'PE':
        end_type_str = 'yes'
    elif ended_type == 'SE':
        end_type_str = 'no'
    else:
        raise ValueError(f'Unrecognized end type {ended_type}')

    out_path = task_context.path_for_volume_file(f'{prefix}_exon_exp.txt')

    # bam_path = task_context.download_s3_object(bam_uri)
    bam_path = task_context.get_result_for_stage('star', 'bam').get_array_result().download_to_volume()
    sorted_bam_path = bam_path.replace('.bam', '.sorted.bam')
    gff_path = task_context.get_asset(f'gencode.v31.annotation.gff')
    # gff_path = task_context.get_asset('gencode.v31.annotation.nochr.gff')
    # gff_path = task_context.get_asset(f'gencode.v31.annotation_SYNGAP_exon31_alt.gff')
    # gff_path = task_context.get_asset(f'gencode.vM25.chr_patch_hapl_scaff.annotation.gff')
    # gff_path = task_context.get_asset(f'rheMac8.ensGene.gff')
    sys.stdout.flush()
    print(f':::INDEXING {prefix}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'samtools sort {bam_path} > {sorted_bam_path}'
        ),
        mappings={
            'bam_path' : bam_path,
            'sorted_bam_path' : sorted_bam_path
        }
    )

    task_context.run_template(
        template=(
            'samtools index {bam_path}'
        ),
        mappings={
            'bam_path' : sorted_bam_path
        }
    )


    sys.stdout.flush()
    print(f':::GETTING COUNTS FOR {prefix}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'python3.6 /home/dexseq_count_original.py '
            '-p {end_type} '
            '-s no '
            '-f bam '
            '-r pos '
            '{gff_path} '
            '{bam_path} '
            '{out_path}'
        ),
        mappings={
            'end_type' : end_type_str,
            'gff_path' : gff_path,
            'bam_path' : sorted_bam_path,
            'out_path' : out_path
        }
    )

    sys.stdout.flush()
    print(f':::FINISHED; UPLOADING {out_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('exon_expression', 'txt', out_path)
    task_context.put_result_for_stage('sorted_bam', 'bam', sorted_bam_path)
    task_context.put_result_for_stage('sorted_bam', 'bai', f'{sorted_bam_path}.bai')

    sys.stdout.flush()

if __name__ == '__main__':
    run_exon_expression()
