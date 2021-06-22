#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from Primer_range import Primer_range
from os.path import basename,dirname
from collections import defaultdict
from Alignment import Alignment
from Fragment import Fragment
from Bio.Seq import Seq
import argparse
import gzip
import sys
import re
import os

parser = argparse.ArgumentParser(description='Script adapted from https://github.com/ItokawaK/Alt_nCov2019_primers, This tool trims portion of primer from aligned reads. '
                                              'Input name sorted sam or output of the bwa mem directry by PIPE.')
parser.add_argument('primer_bed', help='bed file describing primer coordinates')
parser.add_argument('ref', help='fasta file of reference genome')
parser.add_argument('out', help=' output prefix, script will append the amplicon region name')
parser.add_argument('--table',default="",help="output prefix, file storing all pair of amplicon seen and effectif")
parser.add_argument('-s',default=1,help=" 0: only keep properly paired reads, 1: same as 0 and but only keep reads with both primer presents, 2: same as 1 but only keep reads with both primer complete")
parser.add_argument('--log',default='',help="output prefix, stats on filtering")
args = parser.parse_args()


OVERLAP_THRESHOLD = 100
FRAGMENT_MIN_SIZE = 100

# pass args, maybe useless
BED_FILE = args.primer_bed
REF_FILE = args.ref
STRICT = int(args.s)
prefix = args.out


# open 2 handle per amplicon for writing, gziped
primer_range = Primer_range(BED_FILE,REF_FILE)
amplicons = {range[-1].split("_LEFT")[0] for range in primer_range.ranges_fwd}
os.system("mkdir -p %s"%dirname(prefix))
amplicons_to_handles = {amp:[gzip.open("%s_%s_R1.fastq.gz"%(prefix,amp), 'wb', compresslevel=3),gzip.open("%s_%s_R2.fastq.gz"%(prefix,amp), 'wb', compresslevel=3)] for amp in amplicons}


Log = defaultdict(int)
primer_to_cnt=defaultdict(int)
read_to_bucket = defaultdict(lambda :[None,None])

# for index,sam_line in enumerate(open(file)):
for index,sam_line in enumerate(sys.stdin):
    # Skip header
    if sam_line.startswith('@'):
        continue

    alignment = Alignment(sam_line.rstrip())

    # Initialize alignment backet
    current_read = alignment.read_name
    alignment_bucket = read_to_bucket[current_read]

    # Skip non-primary and supplemental alignments
    if alignment.flag & (256 + 2048):
        Log["Not primary or supplementary"]+=1
        continue

    # check for multimapping
    if alignment_bucket[0] and alignment_bucket[1]:
        Log["multimapping reads "]+=1

    if alignment.strand == '+':
        alignment_bucket[0] = alignment
    else:
        alignment_bucket[1] = alignment

    # process if the bucket is full
    if alignment_bucket[0] and alignment_bucket[1]:
        Log[" % properly paired"]+=1
        
        # check that reads overlap, otherwise, it's shit
        overlap = (alignment_bucket[0].ref_start<alignment_bucket[1].ref_end)&(alignment_bucket[1].ref_start-alignment_bucket[0].ref_end<=OVERLAP_THRESHOLD)
        if not overlap:
            Log["Not overlapping"]+=1
            continue
        
        # remove paired reads if fragment length is less than read length
        fragment = Fragment(alignment_bucket[0], alignment_bucket[1])
        if fragment.size<FRAGMENT_MIN_SIZE:
            Log["fragment size < min size"]+=1
            # continue
       
        # get the amplicon region in which reads are contained
        def amp_overlap(f,x):
            start,end,name = x
            Length = float(end-start)
            if start<=f.ref_start<=end:
                return name,(end-f.ref_start)/Length
            elif start<=f.ref_end<=end:
                return name,(f.ref_end-start)/Length
            else:
                return name,0
        try : 
            amp = max([(name,overlap) for name,overlap in map(lambda x:amp_overlap(fragment,x),primer_range.amplicon) if overlap>0.20],key=lambda x:x[1])[0]
        except:
            Log["too small overlap with amplicon"]+=1     
            continue

        # and check if fragment ends are contained in a primer region
        (range_left,amp_left) = primer_range.is_contained(fragment.ref_start, 'left')
        (range_right,amp_right) = primer_range.is_contained(fragment.ref_end, 'right')
       
        if (range_left!=None)&(range_right!=None):
            Log["both primer"]+=1
            # stats on correct primer use : 
            m = primer_range.exact_match(fragment,range_left[-1],range_right[-1])
            if m:
                Log['perfect match with both primer']+=1
            elif STRICT>1:
                continue
        elif (range_left!=None)|(range_right!=None):
            Log["unique primer"]+=1
            if STRICT>0:
                continue
        else:
            Log["no primer"]+=1
            if STRICT>0:
                continue

        primer_to_cnt[amp]+=1
        # Chop ends of the fragment overlapping to primer
        fragment.slice(range_left)
        fragment.slice(range_right)

        # Get fasta string list
        fastq_lines = fragment.get_fastqlines()

        # output reads 
        out_read1 = "%s\n"%"\n".join(fastq_lines[0])
        out_read2 = "%s\n"%"\n".join(fastq_lines[1])
        amplicons_to_handles[amp][0].write(out_read1.encode())
        amplicons_to_handles[amp][1].write(out_read2.encode())

# close thoses handles, otherwise python need to do it himself at the end of the run
# and that cause delay in file appearing on disk, and snakemake fail.
for handle1,handle2 in amplicons_to_handles.values():
    handle1.close()
    handle2.close()

Log["properly paired"]/=float(len(read_to_bucket))
if args.log:
    with open(args.log,"w") as handle:
        sorted_keys = ['Not primary or supplementary',' % properly paired','no primer','unique primer','both primer','wrong primer combination','perfect match with both primer']
        handle.writelines("%s\t%s\n"%(key,Log[key]) for key in sorted_keys)

if args.table:
    with open(args.table,"w") as handle:
        handle.writelines("%s\t%s\n"%(key,values) for key,values in sorted(primer_to_cnt.items(),key=lambda x:x[0]))

