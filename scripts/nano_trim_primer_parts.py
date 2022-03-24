#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from Primer_range import Primer_range
from os.path import basename,dirname
from collections import defaultdict
from Alignment import Alignment
from nano_Fragment import Fragment
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
parser.add_argument('--log',default='',help="output prefix, stats on filtering")
args = parser.parse_args()


FRAGMENT_MIN_SIZE = 100
MIN_RATIO = 0.95
# pass args, maybe useless
BED_FILE = args.primer_bed
REF_FILE = args.ref
prefix = args.out

# open 2 handle per amplicon for writing, gziped
primer_range = Primer_range(BED_FILE,REF_FILE)
amplicons = {range[-1].split("_LEFT")[0] for range in primer_range.ranges_fwd}
os.system("mkdir -p %s"%dirname(prefix))
amplicons_to_handles = {amp:gzip.open("%s_%s.fastq.gz"%(prefix,amp), 'wb', compresslevel=3) for amp in amplicons}


Log = defaultdict(int)
primer_to_cnt=defaultdict(int)
ratios = defaultdict(list)

# debug
# BED_FILE = "/mnt/gpfs/seb/Project/CovMix/primer_schemes/Artic_V4-6/Artic_V4-6_primer.bed"
# REF_FILE = "/mnt/gpfs/seb/Project/CovMix/primer_schemes/Artic_V4-6/Artic_V4-6_reference.fasta"
# file = "/mnt/gpfs/seb/Project/CovMix/nanopore/Airport/samples/A1_propmix2/A1_propmix2_init.sam"

# header,seq = next(sfp(open(REF_FILE)))

# test = []
# handle = open(file)
# sam_line = next(handle)
# for index,sam_line in enumerate(open(file)):
for index,sam_line in enumerate(sys.stdin):
    # Skip header
    if sam_line.startswith('@'):
        continue

    alignment = Alignment(sam_line.rstrip())
    current_read = alignment.read_name

    # Skip non-primary and supplemental alignments
    if alignment.flag & (256 + 2048):
        Log["Not primary or supplementary"]+=1
        continue

    # remove paired reads if fragment length is less than read length
    fragment = Fragment(alignment)
    if fragment.size<FRAGMENT_MIN_SIZE:
        Log["fragment size < min size"]+=1
        # continue
   
    # get the amplicon region in which reads are contained
    def amp_overlap(read,amp):
        start,end,name = amp
        Length = float(end-start)
        # case overlap before amp
        if (read.ref_end<=start)|(end<=read.ref_start):
            return name,0
        return name,(min(end,read.ref_end)-max(start,read.ref_start))/Length
    try : 
        amp = max([(name,overlap) for name,overlap in map(lambda x:amp_overlap(fragment,x),primer_range.amplicon) if overlap>0.20],key=lambda x:x[1])[0]
    except:
        Log["too small overlap with amplicon"]+=1
        continue


    # change the way we trim, do that with just primer definition
    range_left = next((start,end,strand,name) for (start,end,strand,name) in primer_range.ranges_fwd if "alt" not in name if amp in name)
    range_right = next((start,end,strand,name) for (start,end,strand,name) in primer_range.ranges_rev if "alt" not in name if amp in name)
    amp_def = [range_left[1],range_right[0]]


    # Chop ends of the fragment overlapping to primer
    fragment.slice(amp_def)

    # do some stats on indels
    if fragment.ratio!=0:
        Log['both primer']+=1
        if fragment.ratio<MIN_RATIO:
            continue
        ratios[amp].append(fragment.ratio)


    # and check if fragment ends are contained in a primer region
    primer_to_cnt[amp]+=1

    Log["processed succesfully"]+=1


    # Get fasta string list
    fastq_lines = fragment.get_fastqlines()

    # output reads 
    out_read1 = "%s\n"%"\n".join(fastq_lines)
    amplicons_to_handles[amp].write(out_read1.encode())

# close thoses handles, otherwise python need to do it himself at the end of the run
# and that cause delay in file appearing on disk, and snakemake fail.
for handle in amplicons_to_handles.values():
    handle.close()

thresh = lambda x,t:len([el for el in x if el<=t])
ratios_stats= [(amp,str(thresh(rats,0.95)),str(thresh(rats,0.97)),str(thresh(rats,1)),str(thresh(rats,10))) for amp,rats in ratios.items()]

with open("%s_ratios.tsv"%prefix,"w") as handle:
    handle.write("amp\t<=0.95\t<=0.98\t<=1\t<10\n")
    handle.writelines("%s\n"%"\t".join(line) for line in ratios_stats)

if args.log:
    with open(args.log,"w") as handle:
        sorted_keys = ['Not primary or supplementary','processed succesfully','no primer','unique primer','both primer','wrong primer combination','perfect match with both primer']
        handle.writelines("%s\t%s\n"%(key,Log[key]) for key in sorted_keys)

if args.table:
    with open(args.table,"w") as handle:
        handle.writelines("%s\t%s\n"%(key,values) for key,values in sorted(primer_to_cnt.items(),key=lambda x:x[0]))