from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

@pipeline_task
def run_exon_ratios(task_context):
    '''
    Manifest for this stage will be different than main olego manifest
    '''
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix = my_line[0]
    pcd = my_line[1]
    sex = my_line[5]
    bam_uri = my_line[7]
    gene_name = my_line[8]
    regions = my_line[9]
    out_path = task_context.path_for_volume_file(f'{prefix}_{gene_name}.txt')

    sys.stdout.flush()
    print(f':::EXTRACTING REGIONS {regions} FROM {prefix}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'python3.6 /home/get_exon_ratios.py '
            '-s {prefix} '
            '-p {pcd} '
            '-x {sex} '
            '-b {bam_uri} '
            '-r {regions} '
            '-o {out_path}'
        ),
        mappings={
            'prefix' : prefix,
            'pcd': pcd,
            'sex' : sex,
            'bam_uri' : bam_uri,
            'regions' : regions,
            'out_path' : out_path
        }
    )

    sys.stdout.flush()
    print(f':::FINISHED EXTRACTING; UPLOADING {out_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('exon_ratios', 'single_samples', out_path)
    sys.stdout.flush()

if __name__ == '__main__':
    run_exon_ratios()
