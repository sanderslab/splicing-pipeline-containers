'''
A modified version of DEXSeq count (edited by Lindsay Liang; please email lindsay1807@gmail.com for help).

This version of DEXSeq_count will take in an extra 3th positional argument of a list of exons
to 'merge' counts for.  Eg, this script will count the union of reads that align to multiple exons
(the counts for individual exons of a merged exon group will all be 0).

Note that the extra merge functions in this script only work for PE reverse bams. (eg Brainvar)

'''

from __future__ import division
from pprint import pprint
import sys, itertools, optparse, warnings

optParser = optparse.OptionParser(
    
    usage = 'python %prog [options] <flattened_gff_file> <alignment_file> <output_file>',
    
    description=
        'This script COUNTS how many reads in <alignment_file> fall onto each exonic ' +
        'part given in <flattened_gff_file> and outputs a list of COUNTS in ' +
        '<output_file>, for further analysis with the DEXSeq Bioconductor package. ' +
        'Notes: Use dexseq_prepare_annotation.py to produce <flattened_gff_file>. ' + 
        '<alignment_file> may be \'-\' to indicate standard input.',
        
    epilog = 
        'Written by Simon Anders (sanders@fs.tum.de) and Alejandro Reyes (reyes@embl.de), ' +
        'European Molecular Biology Laboratory (EMBL). (c) 2010-2013. Released under the ' +
        ' terms of the GNU General Public License v3. Part of the \'DEXSeq\' package.')
        
optParser.add_option('-p', '--paired', type='choice', dest='paired',
    choices = ('no', 'yes'), default = 'no',
    help = '\'yes\' or \'no\'. Indicates whether the data is paired-end (default: no)')

optParser.add_option('-s', '--stranded', type='choice', dest='stranded',
    choices = ('yes', 'no', 'reverse'), default = 'yes',
    help = '\'yes\', \'no\', or \'reverse\'. Indicates whether the data is ' +
        'from a strand-specific assay (default: yes). ' +
        'Be sure to switch to \'no\' if you use a non strand-specific RNA-Seq library ' +
        'preparation protocol. \'reverse\' inverts strands and is needed for certain ' +
        'protocols, e.g. paired-end with circularization.' )
    
optParser.add_option('-a', '--minaqual', type='int', dest='minaqual',
    default = 10,
    help = 'skip all reads with alignment quality lower than the given ' +
        'minimum value (default: 10)')

optParser.add_option('-f', '--format', type='choice', dest='alignment',
    choices=('sam', 'bam'), default='sam',
    help = '\'sam\' or \'bam\'. Format of <alignment file> (default: sam)')

optParser.add_option('-r', '--order', type='choice', dest='order',
    choices=('pos', 'name'), default='name',
    help = '\'pos\' or \'name\'. Sorting order of <alignment_file> (default: name). Paired-end sequencing ' +
        'data must be sorted either by position or by read name, and the sorting order ' +
        'must be specified. Ignored for single-end data.')
   

def make_merge_sets(merge_set_file):
    '''
    Each line in merge_set_file is two columns - an exon, and the group that it belongs to
    eg:
    A   1
    B   1
    C   2
    D   2

    will merge exons A and B, and C and D separately.

    Each line in merge_set_file is a set of features we want to count as
    single exons; columns are the groups to be merged; must be > 2.
    Each column must follow the format ENSG00000136531.16+ENSG00000236283.5:004
    (gene_id:exon_number)

    '''
    merge_sets = list()
    # For every exon to be merged, map which pair they belong to
    merge_exons = dict()

    exons_to_groups = dict()
    groups_to_exons = dict()

    with open(merge_set_file, 'r') as in_f:
        for line in in_f:
            exon, group = [i.strip() for i in line.split('\t')]
            if group not in groups_to_exons:
                groups_to_exons[group] = [exon]
            else:
                groups_to_exons[group].append(exon)
            if exon not in exons_to_groups:
                exons_to_groups[exon] = group
            else:
                raise ValueError('Every exon can only belong to one group')
    for exon in exons_to_groups:
        group = exons_to_groups[exon]
        exon_group = set(groups_to_exons[group])
        merge_exons[exon] = exon_group
        merge_sets.append(exon_group)

    # with open(merge_set_file, 'r') as in_f:
    #     for line in in_f:
    #         exons = [i.strip() for i in line.split('\t')]
    #         exon_group = set(exons)
    #         merge_sets.append(exon_group)
    #         for exon in exons:
    #             print(merge_exons)
    #             if exon not in merge_exons:
    #                 merge_exons[exon] = exon_group
    #             else:
    #                 # shouldnt overlap
    #                 raise ValueError(f'One exon cannnot belong to more than one group {exon}')
    return merge_sets, merge_exons

#We need this little helper below:
def reverse_strand(s):
    '''
    '''
    if s == '+':
        return '-'
    elif s == '-':
        return '+'
    else:
        raise SystemError('illegal strand')

def update_count_vector(COUNTS, rs):
    '''
    '''
    if(type(rs) == str):
        COUNTS[rs] += 1
    else:
        for f in rs:
            COUNTS[f] += 1
    return COUNTS

def map_read_pair(af, ar):
    '''
    '''
    rs = set()
    if af and ar and not af.aligned and not ar.aligned:
        return '_notaligned'
    if af and ar and not af.aQual < minaqual and ar.aQual < minaqual:
        return '_lowaqual'
    if af and af.aligned and af.aQual >= minaqual and af.iv.chrom in list(features.chrom_vectors.keys()):
        # if the forward read is aligned + passes qual + aligned to one of the chroms in the gff
        for cigop in af.cigar:
            # in the cigar 20M5D9M, 20M, 5D, and 9M are each individual cigops.
            # af.cigar: [< CigarOperation: 101 base(s) matched on ref iv chr2:[164841485,164841586)/+, query iv [0,101) >]
            if cigop.type != 'M':
                # position is not exactly aligned with ref, go to the next op
                continue
            if reverse:
                # if the strandedness argument was 'reverse'
                cigop.ref_iv.strand = reverse_strand(cigop.ref_iv.strand)
            import pdb; pdb.set_trace()
            for iv, s in features[cigop.ref_iv].steps():
                # in the cigar 20M5D9M, 20, 5, and 9 are the ref_iv intervals.
                # for interval, annotation in <ChromVector object, chr2:[165091059,165091154)/+, step>
                # the ChromVector is the 95 bp of the cigar str for this 95bp read; not sure what the steps are...
                # rs is a set of regions from the gff
                rs = rs.union(s)
    if ar and ar.aligned and ar.aQual >= minaqual and ar.iv.chrom in list(features.chrom_vectors.keys()):
        for cigop in ar.cigar:
            if cigop.type != 'M':
                continue
            if not reverse:
                cigop.ref_iv.strand = reverse_strand(cigop.ref_iv.strand)
            for iv, s in features[cigop.ref_iv].steps():
                    # Here's where s isn't being added to rs.
                    rs = rs.union(s)
    set_of_gene_names = set([f.split(':')[0] for f in rs])
    if len(set_of_gene_names) == 0:
        return '_empty'
    elif len(set_of_gene_names) > 1:
        return '_ambiguous'
    else:
        # print(rs)
        return rs

def clean_read_queue(queue, current_position):
    '''
    '''
    clean_queue = dict(queue)
    for i in queue:
        if queue[i].mate_start.pos < current_position:
            warnings.warn('Read '+ i + ' claims to have an aligned mate that could not be found in the same chromosome.')
            del clean_queue[i]
    return clean_queue

def single_end(sam_file, num_reads):
    '''
    Processing reads in single end mode
    '''
    for a in reader(sam_file):
        if not a.aligned:
            print(f'Read {a.read.name} was not aligned')
            COUNTS['_notaligned'] += 1
            continue
        if 'NH' in a.optional_fields and a.optional_field('NH') > 1:
            print(f'Read {a.read.name} had NH > 1')
            continue
        if a.aQual < minaqual:
            print(f'Read {a.read.name} had low qual of {a.aQual}')
            COUNTS['_lowaqual'] += 1
            continue
        rs = set()
        for cigop in a.cigar:
            if cigop.type != 'M':
                print(f'Cigop != M for a {a.read.name}, {cigop}')
                continue
            if reverse:
                cigop.ref_iv.strand = reverse_strand(cigop.ref_iv.strand)
            for iv, s in features[cigop.ref_iv].steps():
                print(f'Added {s} to rs for ar; rs is now {rs}')
                print(f'Cigop iv: {iv}')
                rs = rs.union(s)
        set_of_gene_names = set([f.split(':')[0] for f in rs])
        print(f'Set of gene names for {a.read.name}: {set_of_gene_names}')
        if len(set_of_gene_names) == 0:
            COUNTS['_empty'] += 1
        elif len(set_of_gene_names) > 1:
            COUNTS['_ambiguous'] += 1
        else:
            for f in rs:
                COUNTS[f] += 1
        num_reads += 1
        if num_reads % 100000 == 0:
            sys.stderr.write('%d reads processed.\n' % num_reads)
    return num_reads

def paired_end_pos_sort(sam_file, num_reads, merge_sets, merge_exons):
    '''
    New implementation of paired end processing where bam is sorted by position.
    (Position sorting is the default of samtools sort)
    '''
    mates = dict()
    count = 0

    for a in reader(sam_file):
        if not a.aligned:
            print(f'Read {a.read.name} was not aligned')
            COUNTS['_notaligned'] += 1
            continue
        if 'NH' in a.optional_fields and a.optional_field('NH') > 1:
            print(f'Read {a.read.name} had NH > 1')
            # if this read is aligned to more than one place
            continue
        if a.aQual < minaqual:
            print(f'Read {a.read.name} had low qual of {a.aQual}')
            COUNTS['_lowaqual'] += 1
            continue
        curr_pos = a.iv.start
        curr_chrom = a.iv.chrom
        if a.read.name and a.mate_aligned:
            # if read has a name and aligned
            # print(f'Read {a.read.name} is mate aligned')
            if a.read.name in mates:
                print(f'Already seen {a.read.name}; pair found.')    
                # get previously seen read that's part of this pair
                b = mates[a.read.name]
                # print(a.pe_which, a)
                # print(b.pe_which, b)
                if a.pe_which == 'first' and b.pe_which == 'second':
                    # assigning forward/reverse reads
                    af=b
                    ar=a
                else:
                    af=a
                    ar=b
                print(f'Forward read: {af}')
                print(f'Reverse read: {ar}')
                rs = set()
                for cigop in af.cigar:
                    print(f'Cigop_f: {cigop}, {cigop.type}{cigop.size}')
                    # in the cigar 20M5D9M, 20M, 5D, and 9M are each individual cigops.
                    if cigop.type != 'M':
                        # print(f'Cigop != M for af {af.read.name}, {cigop}')
                        continue
                    for iv, s in features[cigop.ref_iv].steps():
                        # in the cigar 20M5D9M, 20, 5, and 9 are the ref_iv intervals.
                        # iv - interval
                        # s - set of feature returned by steps() for the interval iv
                        # but our features should never have more than one exon in a step
                        # bc none of the gff regions overlaps; might have overlapping
                        # genes in steps
                        print(f'Cigop_f iv: {iv}')
                        rs = rs.union(s)
                        print(f'Added {s} to rs for af; rs is now {rs}')
                for cigop in ar.cigar:
                    print(f'Cigop_r: {cigop}, {cigop.type}{cigop.size}')
                    if cigop.type != 'M':
                        print(f'Cigop != M for ar {ar.read.name}, {cigop}')
                        continue
                    cigop.ref_iv.strand = reverse_strand(cigop.ref_iv.strand)
                    for iv, s in features[cigop.ref_iv].steps():
                        print(f'Cigop_r iv: {iv}')
                        rs = rs.union(s)
                        print(f'Added {s} to rs for ar; rs is now {rs}')
                count += 1
                # if count == 5:
                #     exit()
                # removing a.read.name from mates frees up ram
                del mates[a.read.name]

                # Checking if a read aligned to not exactly one *gene* (not exon)
                set_of_gene_names = set([f.split(':')[0] for f in rs])
                print(f'Set of gene names for {a.read.name}: {set_of_gene_names}')
                if len(set_of_gene_names) == 0:
                    COUNTS['_empty'] += 1
                elif len(set_of_gene_names) > 1:
                    COUNTS['_ambiguous'] += 1
                else:
                    print(f'Processing rs {rs} for read {a.read.name}')
                    # Find if read overlays any of the merging exons (union)
                    seen_merged = set()
                    for f in rs:
                        if f in merge_exons:
                            exon_group = merge_exons[f]
                            exon_nums = [exon.split(':')[1] for exon in exon_group]
                            exon_nums.sort()
                            gene_id = [exon.split(':')[0] for exon in exon_group][0]
                            merge_name = f'{gene_id}:{"+".join(exon_nums)}'
                            if merge_name not in seen_merged:
                                print(f'### MERGED EXON READ {a.read.name}')
                                if merge_name in COUNTS:
                                    COUNTS[merge_name] += 1
                                else:
                                    COUNTS[merge_name] = 1
                                seen_merged.add(merge_name)
                        else:
                            print(f'Adding +1 to COUNTS at feature {f}')
                            COUNTS[f] += 1

                    # Below is the intersection logic
                    # if len(rs) > 1:
                    #     for merge_pair in merge_sets:
                    #         # Check if rs is a subset of any merging pairs
                    #         if merge_pair - rs == set():
                    #             print(f'Found merge pair {merge_pair} of exons in rs')
                    #             # this merge_pair is covered by rs' read
                                
                    #             # generate new name for merged exon
                    #             merge_exon_num = [f.split(':')[1] for f in merge_pair]
                    #             merge_exon_num.sort()

                    #             # all gene_ids in merge_pair should be the same
                    #             gene_id = [f.split(':')[0] for f in merge_pair][0]
                                
                    #             merge_name = f'{gene_id}:{"+".join(merge_exon_num)}'
                    #             if merge_name in COUNTS:
                    #                 COUNTS[merge_name] += 1
                    #             else:
                    #                 COUNTS[merge_name] = 1
                    #             # remove all features counted from rs so
                    #             # they don't get counted as individual features
                    #             rs = rs - merge_pair
                    #             print(f'Removed merge pair {merge_pair} from rs, rs now {rs}')
                    # # After checking for pairs and removing paired features,
                    # # count remaining features
                    # for f in rs:
                    #     print(f'Adding +1 to COUNTS at feature {f}')
                    #     COUNTS[f] += 1
            else:
                # First time encountering this mate pair
                # print(f'First time encountering {a.read.name}')
                if a.mate_start.chrom != a.iv.chrom:
                    print(f'Mates mapped to different chroms! {a.mate_start.chrom}, {a.iv.chrom}')
                    COUNTS['_ambiguous_readpair_position'] += 1
                    continue
                else:
                    # print(f'Adding {a.read.name} to mates.')
                    mates[a.read.name] = a
        else:
            continue
        # book keeping
        num_reads += 1
        if num_reads % 10000 == 0:
            sys.stderr.write('%d reads processed.\n' % num_reads)
    print(f'All unmated reads at the end of bam parsing:')
    pprint(mates)
    return num_reads

if len(sys.argv) == 1:
    optParser.print_help()
    sys.exit(1)

(opts, args) = optParser.parse_args()

if len(args) != 4:
    sys.stderr.write(sys.argv[0] + ': Error: Please provide three arguments.\n')
    sys.stderr.write('  Call with \'-h\' to get usage information.\n')
    sys.exit(1)

try:
    import HTSeq
except ImportError:
    sys.stderr.write('Could not import HTSeq. Please install the HTSeq Python framework\n')    
    sys.stderr.write('available from http://www-huber.embl.de/users/anders/HTSeq\n')    
    sys.exit(1)

gff_file = args[0]
sam_file = args[1]
merge_set_file = args[2]
out_file = args[3]
stranded = opts.stranded == 'yes' or opts.stranded == 'reverse'
reverse = opts.stranded == 'reverse'
is_PE = opts.paired == 'yes'
alignment = opts.alignment
minaqual = opts.minaqual
order = opts.order

if alignment == 'bam':
    try:
        import pysam
    except ImportError:
        sys.stderr.write('Could not import pysam, which is needed to process BAM file (though\n')
        sys.stderr.write('not to process text SAM files). Please install the \'pysam\' library from\n')
        sys.stderr.write('https://code.google.com/p/pysam/\n')    
        sys.exit(1)



if sam_file == '-':
    sam_file = sys.stdin


# Step 1: Read in the GFF file as generated by aggregate_genes.py
# and put everything into a GenomicArrayOfSets

features = HTSeq.GenomicArrayOfSets('auto', stranded=stranded)
print(features)  
for f in  HTSeq.GFF_Reader(gff_file):
    if f.type == 'exonic_part':
        f.name = f.attr['gene_id'] + ':' + f.attr['exonic_part_number']
        features[f.iv] += f.name


# initialise counters
num_reads = 0
COUNTS = {}
COUNTS['_empty'] = 0
COUNTS['_ambiguous'] = 0
COUNTS['_lowaqual'] = 0
COUNTS['_notaligned'] = 0
COUNTS['_ambiguous_readpair_position'] = 0

# put a zero for each feature ID
for iv, s in features.steps():
    for f in s:
        COUNTS[f] = 0


    
if alignment == 'sam':
    reader = HTSeq.SAM_Reader
else:
    reader = HTSeq.BAM_Reader


# Now go through the aligned reads
num_reads = 0

if not is_PE:
    num_reads = single_end(sam_file=sam_file, num_reads=num_reads)
else: # paired-end
    alignments = dict()
    if order == 'name':
        for af, ar in HTSeq.pair_SAM_alignments(reader(sam_file)):
            if af == None or ar == None:
                continue
            if not ar.aligned:
                continue
            if not af.aligned:
                continue
            elif ar.optional_field('NH') > 1 or af.optional_field('NH') > 1:
                continue
            elif af.iv.chrom != ar.iv.chrom:
                COUNTS['_ambiguous_readpair_position'] += 1
                continue
            else:
                rs = map_read_pair(af, ar)
                COUNTS = update_count_vector(COUNTS, rs)
                num_reads += 1
            if num_reads % 100000 == 0:
                sys.stderr.write('%d reads processed.\n' % num_reads)

    else:
        merge_sets, merge_exons = make_merge_sets(merge_set_file=merge_set_file)
        num_reads = paired_end_pos_sort(
            sam_file=sam_file,
            num_reads=num_reads,
            merge_sets=merge_sets,
            merge_exons=merge_exons)
 

 
# Step 3: Write out the results

fout = open(out_file, 'w')
for fn in sorted(COUNTS.keys()):
    fout.write('%s\t%d\n' % (fn, COUNTS[fn]))
fout.close()
