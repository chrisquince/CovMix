#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from os.path import basename,dirname
from collections import defaultdict
from Bio.Seq import Seq
import argparse
import gzip
import sys
import re
import os

parser = argparse.ArgumentParser(description='Script heavily rewroten from https://github.com/ItokawaK/Alt_nCov2019_primers, This tool trims portion of primer from aligned reads. '
                                              'Input name sorted sam or output of the bwa mem directry by PIPE.')
parser.add_argument('sam', help='name sorted sam file')
parser.add_argument('primer_bed', help='bed file describing primer coordinates')
parser.add_argument('out', help=' output prefix, script will append the amplicon region name')
parser.add_argument('--table',default="",help="output prefix, file storing all pair of amplicon seen and effectif")
parser.add_argument('--log',default='',help="output prefix, stats on filtering")
parser.add_argument('--datatype', choices=["p-reads","ont"],default="p-reads",help="Illumina short paired reads or nanopore long reads")
args = parser.parse_args()

OVERLAP_THRESHOLD = 20
FRAGMENT_MIN_SIZE = 50

# pass args, maybe useless
SAM_FILE = args.sam
BED_FILE = args.primer_bed
prefix = args.out
DATATYPE = args.datatype


# ----- various utils functions -----------

def format_reads(read):
    read_name, flag, ref_name, ref_start, mapq, cigar, next_ref_name, next_ref_start, fragment_len, read_seq, read_q,*_ = read

    # SAM has 1-based coordinate
    ref_start = int(ref_start) - 1
    next_ref_start = int(next_ref_start) - 1
    flag = int(flag)

    # get end position on ref
    # in sam, the ref_start ignore any Soft cliped nucleotide. Meaning that it is possible for the reads to start by 100 nuc, having it's ref_start being 100 after seq start.
    # for simplification of pipeline, we shift back the ref start to take into account real read start.
    # so that we can do the simple substraction: reads start on ref - amplicon start on ref and trim the rigth length
    # H is ignored because the corresponding read seq is not outputed
    # I is ignored because, it will always correspond to region in between match not corresponding to any increment along the ref position
    cig_split = [[int(nb),case] for nb,case in re.findall("(\d+)([A-Z])",cigar) if case not in {'I','H'}]
    cig_index = next(i for i,(nb,case) in enumerate(cig_split) if case=="M")
    shift_start = sum(nb for nb,_ in cig_split[:cig_index])
    ref_start -= shift_start

    # however we need to know where on the ref the alignment ends
    span_al = sum(nb for nb,_ in cig_split)
    ref_end = ref_start + span_al

    # get orientation, 0 mean - and  1 mean + 
    if flag & 16 :
        sign = "-"
    else:
        sign = "+"

    return [read_name,ref_start,ref_end,fragment_len,sign,cigar,read_seq,read_q]


def paired_reads_generator(file,Log):
    # exploit the fact the file is sorted by name to output paired reads
    names = defaultdict(int)
    nb_reads = set()
    with open(file) as handle:
        reads = []
        for line in handle:
            if line[0]=="@":
                continue
            sline = line.rstrip().split("\t")
            flag = int(sline[1])
            name = sline[0]
            names[name]+=1
            nb_reads.add(sline[0])
            Log["Total nb of reads"]=len(nb_reads)
            if flag & (256 + 2048):
                Log["Not primary or supplementary"]+=1
                continue
            if flag & 4:
                Log["Unmapped"]+=1
                continue
            # check for multimapping
            if names[name]>2:
                Log["Reads multi-mapped "]+=1
            if reads:
                if sline[0]==reads[0][0]:
                    Log["properly paired"]+=1
                    reads.append(sline)
                    yield sorted([format_reads(r) for r in reads],key=lambda x:x[1]),Log
                    reads = []
                else:
                    yield [format_reads(r) for r in reads],Log
                    reads = [sline]
            else:
                reads.append(sline)


def ont_reads_generator(file,Log):
    nb_reads = set()
    with open(file) as handle:
        for line in handle:
            if line[0]=="@":
                continue
            sline = line.rstrip().split("\t")
            nb_reads.add(sline[0])
            Log["Total nb of reads"]=len(nb_reads)
            flag = int(sline[1])
            if flag & (256 + 2048):
                Log["Not primary or supplementary"]+=1
                continue
            if flag & 4:
                Log["Unmapped"]+=1
                continue
            yield [format_reads(sline)],Log


def paired_reads_overlap(reads,amp_def):
    amp,(_,_,(start,end)) = amp_def
    length = float(end-start)
    r1,r2 = reads
    frag_start = r1[1]
    frag_end = r2[2]
    if start<=frag_start<=end:
        return amp,(end-frag_start)/length
    elif start<=frag_end<=end:
        return amp,(frag_end-start)/length
    elif (frag_start<=start)&(end<=frag_end):
        return amp,1
    else:
        return amp,0


def ont_overlap(read,amp_def):
    amp,(_,_,(start,end)) = amp_def
    length = float(end-start)
    r = reads[0]
    cig = [(int(nb),case) for nb,case in re.findall("(\d+)(\w{1})",r[5]) if case not in {"H","I"}]

    # real start is the position for the first match, yes we undo what we did in read format, but it is easier this way.
    start_lag = sum([nb for nb,case in cig[:next(index for index,(nb,case) in enumerate(cig) if case=="M")]])
    frag_start = r[1]+start_lag

    # real end is the last matching position
    cig_end = cig[::-1]
    end_lag = sum([nb for nb,case in cig_end[:next(index for index,(nb,case) in enumerate(cig_end) if case=="M")]])
    frag_end = r[2]-end_lag 

    if start<=frag_start<=end:
        return amp,(end-frag_start)/length
    elif start<=frag_end<=end:
        return amp,(frag_end-start)/length
    elif (frag_start<=start)&(end<=frag_end):
        return amp,1
    else:
        return amp,0



def trim_primer(read,amp_coord):
    def do_the_trim(read_seq,read_q,side,sign,cigar,trim):
        cig_split = [(int(nb),case) for nb,case in re.findall("(\d+)(\w{1})",cigar)]

        # explanation: depending on the side we trim the cigar must be looked from left or rigth
        # the cigar is related the sequence after it has been orientated to the same dir as the ref
        # the read as stored in the sam file is always orientated as the ref, hence there is no need to consider it's sign here
        dir_cig = (side=="left")*1 - (side=="right")*1

        seq_pos = 0
        for nb,case in cig_split[::dir_cig]:
            if case in {"S","M","X"}:
                trim -=nb
                seq_pos +=nb
                if trim<=0:
                    seq_pos+= trim
                    break
            elif case == "I":
                seq_pos +=nb
            elif case == "D":
                trim-=nb
                if trim<=0:
                    break
            elif case in {"H"}:
                continue
            else:
                print("issue non handled case %s"%case)

        # if dir_cig is -1, we need to reverse it twice to conserve initial seq direction
        trim_read_seq = read_seq[::dir_cig][seq_pos:][::dir_cig]
        trim_read_q = read_q[::dir_cig][seq_pos:][::dir_cig]

        return trim_read_seq,trim_read_q


    [read_name,ref_start,ref_end,fragment_len,sign,cigar,read_seq,read_q] = read
    trim_left = amp_coord[0] - ref_start
    trim_rigth = ref_end - amp_coord[1]
    if trim_left>0:
        read_seq,read_q = do_the_trim(read_seq,read_q,"left",sign,cigar,trim_left)
    if trim_rigth>0:
        read_seq,read_q = do_the_trim(read_seq,read_q,"right",sign,cigar,trim_rigth)

    return read_name,read_seq,read_q,sign


tab = str.maketrans("ACTG", "TGAC")
def reverse_complement_table(seq):
    #https://bioinformatics.stackexchange.com/questions/3583/what-is-the-fastest-way-to-get-the-reverse-complement-of-a-dna-sequence-in-pytho
    return seq.translate(tab)[::-1]

def fastq_reads(read,DATATYPE):
    header,seq,qual,sign = read
    # in sam/bam output, read sequence is reorientated to + by default
    # for ont, if on reverse strand vsearch won't be able to work.
    if DATATYPE!="ont":
        if sign=="-":
            qual = qual[::-1]
            seq = reverse_complement_table(seq)
    return ("@%s\n%s\n+\n%s\n"%(header,seq,qual)).encode()
    


# ----- Real start ot the script -----------

# --- get primer/amplicon definition --- 
# for each amplicon name, 
# - first element is forward primer def
# - second is reverse primer def
# - third is amplicon interval
amp_to_coord = defaultdict(lambda:[[],[],[]])
with open(BED_FILE) as bed_h:
    for line in bed_h:
        ref,start,end,name,_,strand = line.rstrip().split('\t')
        if "alt" in name:
            continue
        amp = name.split("_LEFT")[0].split("_RIGHT")[0]
        if strand == "+":
            amp_to_coord[amp][0] = [int(start), int(end), strand, name]
        if strand == "-":
            amp_to_coord[amp][1] = [int(start), int(end), strand, name]
for amp_name,(fwd,rev,amp) in amp_to_coord.items():
        amp_to_coord[amp_name][2] = [fwd[1],rev[0]]


# iterate over reads/filter/trim and write them

Log = defaultdict(int)
primer_to_cnt=defaultdict(int)

if DATATYPE=="p-reads":
    # generator yielding primary matched paired reads, slighly formated
    read_generator = paired_reads_generator(SAM_FILE,Log)
    # open 2 handle per amplicon for writing, gziped
    amplicons_to_handles = {amp:[gzip.open("%s_%s_R1.fastq.gz"%(prefix,amp), 'wb', compresslevel=3),gzip.open("%s_%s_R2.fastq.gz"%(prefix,amp), 'wb', compresslevel=3)] for amp in amp_to_coord}
    amp_overlap = paired_reads_overlap
else:
    # generator yielding primary matched reads, slighly formated
    read_generator = ont_reads_generator(SAM_FILE,Log)
    # open 1 handle per amplicon for writing, gziped
    amplicons_to_handles = {amp:[gzip.open("%s_%s.fastq.gz"%(prefix,amp), 'wb', compresslevel=3)] for amp in amp_to_coord}
    amp_overlap = ont_overlap


for index,(reads,_) in  enumerate(read_generator):

    # get the amplicon region in which reads are contained
    try : 
        amp = max([(name,overlap) for name,overlap in map(lambda x:amp_overlap(reads,x),amp_to_coord.items()) if overlap>0.20],key=lambda x:x[1])[0]
    except:
        Log["overlap with amplicon < 20%"]+=1     
        continue


    if DATATYPE=="p-reads":
        # check that reads overlap, otherwise, it's shit
    
        # remember reads = [read1,read2]
        # read1 = [name,ref_start,ref_end,frag_len,sign]
        overlap = (reads[0][2]-reads[1][1])>=OVERLAP_THRESHOLD
        if not overlap:
            Log["not overlapping with eachover"]+=1
            continue

        # remove reads if fragment length is less than min size
        frag_len = reads[1][2]-reads[0][1]
        if frag_len<FRAGMENT_MIN_SIZE:
            Log["fragment size < %s "%FRAGMENT_MIN_SIZE]+=1
            continue

    primer_to_cnt[amp]+=1

    # loop so that ont and p-reads can be dealed with in the same fashion
    for index,handle in enumerate(amplicons_to_handles[amp]):
        # Chop ends of the fragment overlapping to primer
        R = trim_primer(reads[index],amp_to_coord[amp][2]) 

        # format reads as fastq
        R = fastq_reads(R,DATATYPE)

        # output reads 
        handle.write(R)

# close thoses handles, otherwise python need to do it himself at the end of the run
# and that cause delay in file appearing on disk, and snakemake fail.
for handles in amplicons_to_handles.values():
    for handle in handles:
        handle.close()

Log["nb valid reads"] = sum(primer_to_cnt.values())
if args.log:
    with open(args.log,"w") as handle:
        sorted_keys = ['Total nb of reads',"nb valid reads",'Not primary or supplementary','properly paired',"overlap with amplicon < 20%","not overlapping with eachover","fragment size < %s "%FRAGMENT_MIN_SIZE,"Unmapped",]
        handle.writelines("%s\t%s\n"%(key,Log[key]) for key in sorted_keys)

if args.table:
    with open(args.table,"w") as handle:
        handle.writelines("%s\t%s\n"%(key,values) for key,values in sorted(primer_to_cnt.items(),key=lambda x:x[0]))

