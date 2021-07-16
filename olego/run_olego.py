from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task, get_data_layer
import os
import sys
import multiprocessing

def __combine_lanes(task_context, prefix, r1_paths, r2_paths):
    '''
    Combine multiple lanes for a sample's R1 or R2
    '''
    combined_r1_path = task_context.path_for_volume_file(f'{prefix}_combined_R1.fastq.gz')
    combined_r2_path = task_context.path_for_volume_file(f'{prefix}_combined_R2.fastq.gz')
    
    sys.stdout.flush()
    print('Running cat on R1')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'cat {paths} > {out}'
        ),
        mappings={
            'paths' : ' '.join(r1_paths),
            'out' : combined_r1_path
        }
    )
    sys.stdout.flush()
    print('Running cat on R2')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'cat {paths} > {out}'
        ),
        mappings={
            'paths' : ' '.join(r2_paths),
            'out' : combined_r2_path
        }
    )
    return combined_r1_path, combined_r2_path

def __run_alignment(r_path, sam_path, ref_path):
    '''
    For multiprocessing
    Takes 103mins on 8 threads
    '''
    task_context = get_data_layer()
    sys.stdout.flush()
    print(f'::RUNNING OLEGO ON {r_path}::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'olego '
            '-o {sam_path} '
            '-t 18 '
            '{ref_path} '
            '{r_path}'
            
        ),
        mappings={
            'r_path': r_path,
            'ref_path': ref_path,
            'sam_path': sam_path
        }
    )

    # task_context.put_result_for_stage('olego', f'sam', sam_path)

@pipeline_task
def olego(task_context):
    '''
    Main method for running olego task
    '''
    ref_suffixes = ['', '.pac', '.ann', '.amb', '.rpac', '.bwt', '.rbwt', '.sa', '.rsa']
    ref_paths = list()

    sys.stdout.flush()
    print('::DOWNLOADING REFERENCES::')
    sys.stdout.flush()
    for suff in ref_suffixes:
        ref_paths.append(task_context.get_asset(f'GRCm38.p6.genome.fa{suff}'))
    ref_path = ref_paths[0]

    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix = my_line[0]
    r_uris = my_line[1:]

    r1_paths = list()
    r2_paths = list()
    for uri in r_uris:
        sys.stdout.flush()
        print('::DOWNLOADING FASTQ::')
        sys.stdout.flush()
        if 'R1' in uri:
            r1_paths.append(task_context.download_s3_object(uri))
        elif 'R2' in uri:
            r2_paths.append(task_context.download_s3_object(uri))
        else:
            raise ValueError(f'Could not determine R1/R2 of {uri}')

    r1_path, r2_path = __combine_lanes(task_context, prefix, r1_paths, r2_paths)
    r1_sam = task_context.path_for_volume_file(f'{prefix}_r1.sam')
    r2_sam = task_context.path_for_volume_file(f'{prefix}_r2.sam')
    combined_sam = task_context.path_for_volume_file(f'{prefix}.sam')

    with multiprocessing.Manager() as manager:
        pool = multiprocessing.Pool(processes=2)
        runner_args = [(r1_path, r1_sam, ref_path), (r2_path, r2_sam, ref_path)]
        pool.starmap(__run_alignment, runner_args)
        sys.stdout.flush()
        print('::FINISHED OLEGO. COMBINING SAMS::')
        sys.stdout.flush()

        task_context.run_template(
          template=(
              'perl /home/mergePEsam_fix.pl '
              '-v '
              '{r1_sam} '
              '{r2_sam} '
              '{combined_sam}'),
          mappings={
              'r1_sam': r1_sam,
              'r2_sam': r2_sam,
              'combined_sam': combined_sam
          }
        )

    task_context.put_result_for_stage('olego', f'combined_sam', combined_sam)

if __name__ == '__main__':
    olego()
