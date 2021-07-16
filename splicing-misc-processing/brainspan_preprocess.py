from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

@pipeline_task
def untar_gzip(task_context):
    '''
    Untar and gzip brainspan fastq
    '''
    # tester(task_context)
    # return
    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    # prefix, bam_uri, ended_type = my_line
    prefix  = my_line[0]
    read_uri = my_line[1]
    ended_type = my_line[-1] # SE vs PE

    if ended_type == 'PE':
        end_type_str = 'yes'
    elif ended_type == 'SE':
        end_type_str = 'no'
    else:
        raise ValueError(f'Unrecognized end type {ended_type}')

    out_path = task_context.path_for_volume_file(f'{prefix}.fastq.gz')

    tar_path = task_context.download_s3_object(read_uri)
    untar_path = tar_path.strip('.tbz')
    
    sys.stdout.flush()
    print(f':::UNTARRING {tar_path}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'tar xjf {tar_path} --directory {parent}'
        ),
        mappings={
            'tar_path' : tar_path,
            'parent' : Path(tar_path).parent
        }
    )

    sys.stdout.flush()
    print(f':::GZIPPING {untar_path}:::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            'gzip -c {untar_path} > {out_path}'
        ),
        mappings={
            'untar_path' : untar_path,
            'out_path' : out_path
        }
    )


    sys.stdout.flush()
    print(f':::FINISHED; UPLOADING {out_path}:::')
    sys.stdout.flush()
    task_context.put_result_for_stage('preprocess', 'gz', out_path)

    sys.stdout.flush()

if __name__ == '__main__':
    untar_gzip()
