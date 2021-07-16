from argparse import ArgumentParser
from pathlib import Path
from rkstr8.clients.container import pipeline_task
import csv
import os
import sys

@pipeline_task
def merge_exon_ratios(task_context):
    '''
    Merge exon count files for each gene run
    '''
    counts_objs = task_context.get_result_for_stage('exon_ratios', 'single_samples')
    counts_paths = [count.download_to_volume() for count in counts_objs]

    gene_to_regions_file = task_context.get_asset('brainvar_regions.tsv')

    # Make dict key'd by gene, value'd by regions of the gene
    sys.stdout.flush()
    print('::MAKING REGION DICT::')
    sys.stdout.flush()
    regions_dict = dict()
    with open(gene_to_regions_file, 'r') as regions_f:
        for line in regions_f:
            gene, region = [i.strip() for i in line.split('\t')]
            regions_dict[gene] = region

    # Make dict key'd by gene, value list of paths
    sys.stdout.flush()
    print('::MAKING PATH DICT::')
    sys.stdout.flush()
    path_dict = dict()
    for count in counts_paths:
        count_path = Path(count)
        count_filename = count_path.name
        gene = count_filename.split('_')[-1].strip('.txt')
        if gene not in path_dict:
            path_dict[gene] = [count]
        else:
            path_dict[gene].append(count)

    for gene in path_dict:
        print(f'::WRITING FILE FOR GENE {gene}::')
        sys.stdout.flush()
        out_path = task_context.path_for_volume_file(f'{gene}_readCounts.txt')
        # Write header
        with open(out_path, 'w') as out_f:
            w = csv.writer(out_f, delimiter='\t')
            regions = regions_dict[gene].split(',')
            w.writerow(['#sample', 'age_pcd'] + regions)
        task_context.run_template(
            template=('cat {paths} >> {out_path}'),
            mappings={
                'paths' : ' '.join(path_dict[gene]),
                'out_path' : out_path
            })
        task_context.put_result_for_stage('final_read_counts', gene, out_path)


if __name__ == '__main__':
    merge_exon_ratios()
