from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import os
import sys

def __dl_unpack_refs(task_context, tar_gz, unpack_dir):
    '''
    Download salmon index tarball and untar
    '''
    print('::UNPACKING REFS::')
    sys.stdout.flush()
    ref_tar_gz = task_context.get_asset(tar_gz)
    asset_dir = Path(ref_tar_gz).parent

    task_context.run_template(
        template='tar xzf {package} --directory {parent}',
        mappings={
            'package': ref_tar_gz,
            'parent': asset_dir
        }
    )
    return asset_dir / unpack_dir 

@pipeline_task
def salmon(task_context):
    prefix, r1_uri, r2_uri = tuple(task_context.get_my_input_line().strip().split('\t'))

    print('::DOWNLOADING FASTQs::')
    sys.stdout.flush()

    r1_path = task_context.download_s3_object(r1_uri)
    r2_path = task_context.download_s3_object(r2_uri)

    print('::DOWNLOADING AND UNPACKING REFS::')
    sys.stdout.flush()

    ref_tar = 'gencode.v31.all_transcripts.tar.gz'
    ref_dir = 'gencode.v31.all_transcripts'
    ref_idx_path = __dl_unpack_refs(task_context, ref_tar, ref_dir)

    out_dir = task_context.path_for_volume_file(prefix)
    out_files = {
        'log' :'logs/salmon_quant.log',
        'quant' : 'quant.sf',
        'cmd_info' : 'cmd_info.json',
        'lib_format_counts' : 'lib_format_counts.json',
        'observed_5p_bias' : 'aux_info/observed_bias.gz',
        'observed_3p_bias' : 'aux_info/observed_bias_3p.gz',
        'expected_5p_bias' : 'aux_info/expected_bias.gz',
        'observed_gc_bias' : 'aux_info/obs_gc.gz',
        'expected_gc_biac' : 'aux_info/exp_gc.gz',
        'eq_classes' : 'aux_info/eq_classes.txt.gz',
        'meta_info' : 'aux_info/meta_info.json',
        'fld' : 'aux_info/fld.gz',
        'ambig_info' : 'aux_info/ambig_info.tsv',
        'flenDist' : 'libParams/flenDist.txt'
        }

    print('::RUNNING SALMON::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            'salmon quant -i {ref_idx_path} '
            '-l A '
            '-1 {r1_path} '
            '-2 {r2_path} '
            '-o {out_dir} '
            '--validateMappings '
            '--gcBias '
            '--dumpEq'
        ),
        mappings={
            'ref_idx_path' : ref_idx_path,
            'r1_path' : r1_path,
            'r2_path' : r2_path,
            'out_dir' : out_dir
        }
    )

    print('::UPLOADING RESULTS::')
    sys.stdout.flush()
    for file_type in out_files:
        out_file = f'{out_dir}/{out_files[file_type]}'
        task_context.put_result_for_stage('salmon', file_type, out_file)


if __name__ == '__main__':
    salmon()
