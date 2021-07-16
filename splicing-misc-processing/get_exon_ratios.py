# Description: Find the ratios of different exons in BrainVar data over time
# Usage: python3 olego_brainvar.py chr2:165308703-165308794 chr2:165309165-165309256 chr2:165309352-165309443 chr2:165310323-165310414 
# Author: Stephan Sanders
from argparse import ArgumentParser
import csv
import os
import re
import subprocess
import sys

def printv(x):
    '''
    Quick way of printing variable name and variable, useful when debugging
    '''
    import inspect
    frame = inspect.currentframe().f_back
    s = inspect.getframeinfo(frame).code_context[0]
    r = re.search(r"\((.*)\)", s).group(1)
    print("{} = {}".format(r,x))

def get_ratios(samp_name, pcd, sex, bam_uri, region_list, out_path):
    '''
    Writes a file reporting 
    '''

	# Get the unique reads for the region
    bam_dict = dict()
    region_count = list()
    region_starts = list()
    region_ends = list()

    region_count.append(samp_name)
    region_count.append(pcd)

    for region in region_list:
        region_pos = region.split(':')[1]
        region_start, region_end = [int(i) for i in region_pos.split('-')]
        region_starts.append(region_start)
        region_ends.append(region_end)
        region_count.append(0)
        count = subprocess.run(['samtools', 'view', bam_uri, region], stdout=subprocess.PIPE)
        bam_data = count.stdout.decode('utf-8').split('\n')
        for bam_line in bam_data:
            bam_dict[bam_line] = 5

	# For each read work out the regions with data
    for bam_line in bam_dict:
        bam_fields = bam_line.strip().split('\t')
        if len(bam_fields) > 6:
            start_pos = int(bam_fields[3])
            cigar = bam_fields[5]
            # printv(start_pos)
            # printv(cigar)
            end = 0
            exon_starts = [start_pos]
            exon_ends = list()

            # Find the start and ends of the exons
            while end == 0:
                z = re.search(r'^(\d+)(M|N)', cigar)
                if z:
                    num = int(z.group(1))
                    let = z.group(2)
                    if let == 'M':
                        # match
                        exon_ends.append(exon_starts[-1] + num)
                        cigar = re.sub(r'^(\d+)(M|N)', '', cigar)
                    elif let == 'N':
                        # skip
                        exon_starts.append(exon_ends[-1] + num)
                        cigar = re.sub(r'^(\d+)(M|N)', '', cigar)
                elif cigar == '':
                    end = 1
                else:
                    print(f'ERROR: unexpexted CIGAR string: {cigar}')
                    exon_starts = list()
                    exon_ends = list()
                    end = 1

            # printv(exon_starts)
            # printv(exon_ends)

            # For each region assess overlap with each exon
            for i in range(len(region_list)):
                start_reg = region_starts[i]
                end_reg = region_ends[i]

                # for each exon
                for j in range(len(exon_starts)):
                    start_exon = exon_starts[j]
                    end_exon = exon_ends[j]
                    if start_reg <= end_exon and end_reg >= start_exon:
                        region_count[i+2] = region_count[i+2] + 1


    with open(out_path, 'w') as out_f:
        w = csv.writer(out_f, delimiter='\t')
        print(f'REGION COUNT: {region_count}')
        w.writerow(region_count)


def parse_args():
    parser = ArgumentParser()
    parser.add_argument('-s', '--sample', type=str, required=True)
    parser.add_argument('-p', '--pcd', type=str, required=True)
    parser.add_argument('-x', '--sex', type=str, required=True)
    parser.add_argument('-b', '--bam', type=str, required=True)
    parser.add_argument('-r', '--regions', type=str, required=True)
    parser.add_argument('-o', '--out', type=str, required=True)

    args = parser.parse_args()
    region_list = args.regions.split(',')
    return args.sample, args.pcd, args.sex, args.bam, region_list, args.out

if __name__ == '__main__':
    samp_name, pcd, sex, bam_uri, region_list, out_path = parse_args()
    get_ratios(samp_name, pcd, sex, bam_uri, region_list, out_path)
