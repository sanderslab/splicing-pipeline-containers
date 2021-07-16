from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import csv
import os
import sys


@pipeline_task
def make_junc(task_context):
    '''
    Leafcutter Step 1:
    For each sample, convert olego .sam file to .junc file
    '''

    my_line = tuple(task_context.get_my_input_line().strip().split('\t'))
    prefix = my_line[0]

    sys.stdout.flush()
    print('::DOWNLOADING SAM::')
    sys.stdout.flush()
    sam_path = task_context.get_result_for_stage('olego', 'combined_sam').get_array_result().download_to_volume()
    junc_path = task_context.path_for_volume_file(f'{prefix}.junc')

    sys.stdout.flush()
    print('::MAKING JUNC FILE::')
    sys.stdout.flush()
    task_context.run_template(
        template=(
            '/home/leafcutter/scripts/bam2junc.sh '
            '{sam_path} '
            '{junc_path} '
        ),
        mappings={
            'sam_path': sam_path,
            'junc_path': junc_path
        }
    )
    sys.stdout.flush()
    print('::FINISHED MAKING JUNC. UPLOADING.::')
    sys.stdout.flush()

    task_context.put_result_for_stage('leafcutter_junc', 'junc_file', junc_path)

def __make_junc_path_file(task_context, junc_paths):
    '''
    Leafcutter Step 2:
    Make file listing .junc file paths Step 1

    returns path to junc_file_path file
    '''
    all_juncs_path = task_context.path_for_volume_file('all_juncs.txt')

    sys.stdout.flush()
    print('::MAKING JUNC PATH::')
    sys.stdout.flush()
    with open(all_juncs_path, 'w') as out_f:
        w = csv.writer(out_f, delimiter='\t')
        for path in junc_paths:
            print(f'::JUNC PATH {path}::')
            sys.stdout.flush()
            w.writerow([path])
    
    return all_juncs_path

def __cluster_introns(task_context, all_junc_path):
    '''
    python ../clustering/leafcutter_cluster.py -j test_juncfiles.txt -m 50 -o testYRIvsEU -l 500000
    '''
    sys.stdout.flush()
    print('::CLUSTERING INTRONS::')
    sys.stdout.flush()

    out_gz_prefix = 'cohort_counts'
    out_perind_gz_path = task_context.path_for_volume_file(f'{out_gz_prefix}_perind_numers.counts.gz')
    out_gz_path = task_context.path_for_volume_file(f'{out_gz_prefix}_perind.counts.gz')
    
    parent = Path(out_gz_path).parent

    task_context.run_template(
        template=(
            'python /home/leafcutter/clustering/leafcutter_cluster.py '
            '-j {all_junc_path} '
            '-m 50 ' # 50 split reads supporting each cluster
            '-r {out_dir} '
            '-o {out_gz_prefix} '
            '-l 500000' # allows introns of up to 500kb
        ),
        mappings={
            'all_junc_path': all_junc_path,
            'out_dir' : parent,
            'out_gz_prefix': out_gz_prefix
        }
    )
    sys.stdout.flush()
    print(f'DIR CONTENTS OF {parent} : {os.listdir(parent)}')
    print('::FINISHED INTRON CLUSTERING::')
    task_context.put_result_for_stage('leafcutter_main', 'intron_cluster_tar', out_gz_path)
    task_context.put_result_for_stage('leafcutter_main', 'intron_cluster_tar', out_perind_gz_path)
    sys.stdout.flush()
    return out_gz_path, out_perind_gz_path

def __differential_cluster_analysis(task_context, intron_cluster_tar, gtf_path, group_file_path):
    '''
    ../scripts/leafcutter_ds.R --num_threads 4 --exon_file=../leafcutter/data/gencode19_exons.txt.gz ../example_data/testYRIvsEU_perind_numers.counts.gz example_geuvadis/groups_file.txt
    
    returns {prefix}_cluster_significance.txt and 
    {prefix}_effect_sizes.txt paths
    '''
    out_prefix = task_context.path_for_volume_file('cohort_leafcutter_ds')

    sys.stdout.flush()
    print('::RUNNING DIFFERENTIAL ANALYSIS::')
    sys.stdout.flush()

    task_context.run_template(
        template=(
            '/home/leafcutter/scripts/leafcutter_ds.R '
            '--num_threads {num_threads} '
            '--exon_file={gencode_gtf} '
            '-o {out_prefix} '
            '{intron_cluster_tar} '
            '{group_file}'
        ),
        mappings={
            'num_threads': '8',
            'gencode_gtf': gtf_path,
            'intron_cluster_tar' : intron_cluster_tar,
            'group_file' : group_file_path,
            'out_prefix' : out_prefix
        }
    )


    sig_out = f'{out_prefix}_cluster_significance.txt'
    eff_out = f'{out_prefix}_effect_sizes.txt'

    sys.stdout.flush()
    print('::FINISHED DIFFERENTIAL CLUSTER ANALYSIS. UPLOADING::')
    sys.stdout.flush()

    return sig_out, eff_out

def __plot_leafcutter(task_context, gtf_path, intron_cluster_tar, group_file_path, cluster_significance_path):
    '''
    ./scripts/ds_plots.R -e ../leafcutter/data/gencode19_exons.txt.gz ../example_data/testYRIvsEU_perind_numers.counts.gz example_geuvadis/groups_file.txt leafcutter_ds_cluster_significance.txt -f 0.05
    
    returns path to plot pdf
    '''

    out_path = task_context.path_for_volume_file('cohort_ds_plots.pdf')
    task_context.run_template(
        template=(
            '/home/leafcutter/scripts/ds_plots.R '
            '-e {gencode_gtf} '
            '-o {out_path} '
            '{intron_cluster_tar} '
            '{group_file} '
            '{cluster_significance_path} '
            '-f 0.05'
        ),
        mappings={
            'num_threads': '12',
            'gencode_gtf': gtf_path,
            'intron_cluster_tar' : intron_cluster_tar,
            'group_file' : group_file_path,
            'cluster_significance_path' : cluster_significance_path,
            'out_path' : out_path
        }
    )
    return out_path

@pipeline_task
def leafcutter(task_context):
    '''
    Main leafcutter method
    '''
    sys.stdout.flush()
    print('::DOWNLOADING JUNCS::')
    sys.stdout.flush()

    juncs = task_context.get_result_for_stage('leafcutter_junc', 'junc_file')
    junc_paths = [junc.download_to_volume() for junc in juncs]

    sys.stdout.flush()
    print('::DOWNLOADING GTF AND GROUP FILE::')
    sys.stdout.flush()

    group_file_path = task_context.get_asset('Kriegstein_case_control_groups.txt')
    gtf_path = task_context.get_asset('gencode.v31.exons.txt')

    junc_path_file = __make_junc_path_file(task_context, junc_paths)
    cluster_tar_for_fasqtl, intron_cluster_tar = __cluster_introns(task_context, junc_path_file)
    cluster_sig, eff_size = __differential_cluster_analysis(task_context, intron_cluster_tar, gtf_path, group_file_path)
    plot_pdf = __plot_leafcutter(task_context, gtf_path, intron_cluster_tar, group_file_path, cluster_sig)

    sys.stdout.flush()
    print('::UPLOADING RESULTS::')
    sys.stdout.flush()
    task_context.put_result_for_stage('leafcutter_main', 'intron_cluster_tar', intron_cluster_tar)
    task_context.put_result_for_stage('leafcutter_main', 'intron_cluster_tar', cluster_tar_for_fasqtl)
    task_context.put_result_for_stage('leafcutter_main', 'cluster_significance', cluster_sig)
    task_context.put_result_for_stage('leafcutter_main', 'cluster_effect_size', eff_size)
    task_context.put_result_for_stage('leafcutter_main', 'plot', plot_pdf)
    sys.stdout.flush()
    print('::DONE UPLOADING::')
    sys.stdout.flush()


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-p', '--program', type=str, required=True)
    args = parser.parse_args()
    return args.program

if __name__ == '__main__':
    PROGRAM = parse_args()
    if PROGRAM == 'make_junc':
        make_junc()
    elif PROGRAM == 'main':
        leafcutter()
