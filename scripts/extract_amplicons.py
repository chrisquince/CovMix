#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict,Counter
from multiprocessing import Pool
from os.path import dirname
import numpy as np
import argparse
import glob
import time
import os
import re


parser = argparse.ArgumentParser(description='Extract amplicons from .msa of database')
parser.add_argument('primer_bed', help='bed file describing primer coordinates, beware this script assume some naming conventions, (_LEFT,_RIGHT)')
parser.add_argument('msa', help='fasta file of genome databases')
parser.add_argument('amp', help='prefix for amplicon database sequences')
parser.add_argument('-r', help='fasta file with reference fasta used for primer definition',default="MN908947.3")
parser.add_argument('-s', help='output some stats on database: nb unique amplicons/max nb of shared amplicons ....')
parser.add_argument('-tmin', help='minimum ref position to call variant',default="300")
parser.add_argument('-tmax', help='minimum ref position to call variant',default="29700")
args = parser.parse_args()

MSA_FILE = args.msa
REFNAME = args.r
trim = [float(args.tmin),float(args.tmax)]
BED_FILE = args.primer_bed
AMP_OUT = args.amp
STAT_OUT = args.s

def get_pos_index_map(refseq):
    # here refseq is the reference sequence from a .msa file
    index_to_pos = {}
    pos_to_indexes = defaultdict(list)

    # the add1 is necessary for consistency with variants build under:
    # when looking at variant we concat all base observed until the pos for which ref is != "-". 
    # we then use this pos for the concat bases.
    # here we use the variable add1 for the same purpose
    pos=-1
    add1=0
    for index,nuc in enumerate(refseq):
        if (nuc=="-")&(add1==0):
            add1 = 1
        if nuc!="-":
            pos+=1
            add1=0
        index_to_pos[index]=pos+add1
        pos_to_indexes[pos+add1].append(index)
    return index_to_pos,pos_to_indexes

def get_variant(name_to_seq,refseq,index_to_pos,trim):
    # iterate over the seqs, compare to ref and concat insertions, does not care to concat deletions.
    variants = []
    for name,seq in name_to_seq.items():
        temp=""
        for index,nuc in enumerate(refseq):
            pos = index_to_pos[index]
            if trim[0]<=pos<=trim[1]:
                var_nuc = seq[index]
                if temp:
                    temp+=var_nuc
                    if nuc!="-":
                        if var_nuc==nuc:
                            variants.append([name,index_to_pos[index],temp])
                            temp=""
                else:
                    if nuc!="-":
                        if nuc!=var_nuc:
                            variants.append([name,index_to_pos[index],var_nuc])
                    else:
                        if nuc!=var_nuc:
                            temp+=var_nuc
    return variants


def get_ampseq(amp_to_coord,seq):
    ampseq = {}
    degenerate = {}
    for amp,(start,end) in amp_to_coord.items():
        _ampseq = seq[start:end].upper()
        ampseq[amp]=_ampseq
        if set(_ampseq)-{"A","T","G","C"} != set():
            degenerate[amp]=""
    return ampseq,degenerate


def get_amp_coord(BED_FILE,pos_to_indexes):
    primer_to_coord = defaultdict(lambda:[[],[]])
    amp_to_coord = {}
    with open(BED_FILE) as bed_h:
        for line in bed_h:
            ref,start,end,name,_,strand = line.rstrip().split('\t')
            if "alt" in name:
                continue
            amp = name.split("_LEFT")[0].split("_RIGHT")[0]
            if strand == "+":
                primer_to_coord[amp][0] = [int(start), int(end), strand, name]
            if strand == "-":
                primer_to_coord[amp][1] = [int(start), int(end), strand, name]
    for amp_name,(fwd,rev) in primer_to_coord.items():
            amp_to_coord[amp_name] = [min(pos_to_indexes[fwd[1]]),max(pos_to_indexes[rev[0]])]
    return amp_to_coord


# ------------ real start ------------
name_to_seq = {name:seq for name,seq in sfp(open(MSA_FILE))}
REFSEQ = name_to_seq[REFNAME]

# get column to position mapping
index_to_pos,pos_to_indexes = get_pos_index_map(REFSEQ)

# get amplicon coord
amp_to_coord = get_amp_coord(BED_FILE,pos_to_indexes)

# get primer sequences
header_to_ampseq = {}
header_to_degenerate = {}
for header,seq in name_to_seq.items():
    header = header.split(" ")[0]
    header_to_ampseq[header],header_to_degenerate[header] = get_ampseq(amp_to_coord,seq)

# header_to_degenerate = {header:ampseq for header,ampseq in header_to_degenerate.items() if ampseq!={}}
# nb_degenerate = sum([len(val) for val in header_to_degenerate.values()])

# dereplicate amplicons
amp_to_seq_header = defaultdict(lambda:defaultdict(list))
for header,amp_seq in header_to_ampseq.items():
    for amp,seq in amp_seq.items():
        amp_to_seq_header[amp][seq].append(header)

# renumber/rename amplicons 
ampseq_to_name = {}
ampheader_name = {}
new_name_to_ref = defaultdict(list)
for amp,seq_header in amp_to_seq_header.items():
    for index,(seq,headers) in enumerate(seq_header.items()):
        new_name = '%s_var%s'%(amp,index)
        ampseq_to_name[(amp,seq)] = new_name
        for header in headers:
            ampheader_name[(amp,header)] = new_name 
            new_name_to_ref[new_name].append(header)

# output seqs
os.system("mkdir -p %s"%AMP_OUT)

for amp,seq_header in amp_to_seq_header.items():
    with open("%s/%s.fa"%(AMP_OUT,amp),"w") as handle:
        handle.writelines(">%s\n%s\n"%(ampseq_to_name[(amp,seq)],seq) for (seq,_) in seq_header.items())

# output name mapping
with open("%s/amp_name_mapping.tsv"%(AMP_OUT),"w") as handle:
    handle.writelines("%s\t%s\t%s\n"%(amp,name,header) for (amp,header),name in ampheader_name.items())



# ----------- BONUS look at db ambiguity --------------------

# tabulate unique combination of amplicons and link it to corresponding refs
amp_index = {amp:index for index,amp in enumerate(sorted(amp_to_seq_header.keys()))}
ref_to_combi = defaultdict(lambda:[None for _ in range(len(amp_index))])
for (amp,header),name in ampheader_name.items():
    ref_to_combi[header][amp_index[amp]]=name
combi_to_refs = defaultdict(list)
for ref,combi in ref_to_combi.items():
    combi_to_refs[tuple(combi)].append(ref)


# output combinations
with open("%s/amp_combination_to_ref.tsv"%(AMP_OUT),"w") as handle:
    handle.writelines("%s\t%s\n"%("|".join(combi),"|custom_separator|".join(refs)) for combi,refs in combi_to_refs.items())
    


# nb of unique amplicon per ref,
header_to_uniq = defaultdict(int)
for val in new_name_to_ref.values():
    if len(val)==1:
        header_to_uniq[val[0]]+=1

# nb of amplicons
header_to_nb_amp = {header:len(ampseq) for header,ampseq in header_to_ampseq.items()}

# max nb of amplicons shared with another ref
header_to_shared = {}
for header,ampseq in header_to_ampseq.items():
    header_to_shared[header] = max(Counter([h for amp in ampseq for h in new_name_to_ref[ampheader_name[(amp,header)]] if h!=header]).items(),key=lambda x:x[1])

#output stats
with open("%s/db_stat.tsv"%STAT_OUT,"w") as handle:
    handle.write("name\tnb_amp\tnb_unique\tmax_shared_amp_nb\tmax_shared_ref\n")
    handle.writelines("%s\t%s\t%s\t%s\t%s\n"%(header,header_to_nb_amp[header],header_to_uniq[header],header_to_shared[header][1],header_to_shared[header][0]) for header in header_to_nb_amp)


# --------------- BONUS output Variant -----------------
variants = get_variant(name_to_seq,REFSEQ,index_to_pos,trim)
with open("%s/db_snv.tsv"%dirname(AMP_OUT),"w") as handle:
    handle.writelines("%s\t%s\t%s\n"%(name,str(pos),nuc) for name,pos,nuc in variants)



