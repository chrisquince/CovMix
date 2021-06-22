#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict
from os.path import basename
import argparse
import glob
import re

# REF = "/mnt/gpfs/seb/Database/covid/nCoV-2019.reference.fasta"
# PRIMER = "/mnt/gpfs/seb/Database/covid/nCoV-2019.primer.bed"
# REP = "/mnt/gpfs/seb/Database/covid/GISAID/reps_and_ref.fa"


parser = argparse.ArgumentParser(description='from multiple alignement fasta generate')
parser.add_argument('amp', help='folder with msa of amplicons')
parser.add_argument('map', help='mapping file for amp variant name to refferences')
parser.add_argument('out', help='name of output file')
parser.add_argument('-r',default="MN908947.3",help='name of reference sequence, should be included in the .msa')
args = parser.parse_args()

# pass args, maybe useless
AMP = args.amp
MAP = args.map
OUT = args.out
REF = args.r

# COVID_DB = "/home/sebr/seb/Database/covid"
# COVID_DB_REP = "%s/GISAID/reps_and_ref.fa"%COVID_DB
# AMP,MAP,OUT,REF = "%s/GISAID/amp_seqs"%COVID_DB,"%s/GISAID/amp_seqs/amp_name_mapping.tsv"%COVID_DB,COVID_DB_REP.replace(".fa","_variants.fa"),"MN908947.3"

# amp files
files = glob.glob("%s/*.msa"%AMP)

# name to ref
amp_name_to_headers = defaultdict(list)
for line in open(MAP):
    amp,name,header = line.rstrip().split("\t")
    amp_name_to_headers[(amp,name)].append(header)

# get amp to seq to header for 
header_to_ampseq = defaultdict(lambda:defaultdict(list))
for file in files:
    amp = basename(file.replace(".fa",""))
    for name,seq in sfp(open(file)):
        for header in amp_name_to_headers[(amp.replace(".msa",""),name)]:
            header_to_ampseq[header][amp]=seq

# get variant position for each amplicons
amp_to_var={}
ref_to_var = defaultdict(lambda:defaultdict(list))
for amp,ref in header_to_ampseq[REF].items():
    header_to_variant={}
    for header,ampseq in header_to_ampseq.items():
        if header==ref:
            continue
        seq = ampseq[amp]
        # I don't want to deal with N or stupid Y/K/... 
        if len(set(seq))>4:
            continue
        tmp = []
        for index,nuc in enumerate(seq):
            if nuc!=ref[index]:
                tmp.append((index,nuc))
        header_to_variant[header] = tmp
        if tmp:
            ref_to_var[header][amp] = tmp
    variant_to_refs = defaultdict(list)
    for header,var in header_to_variant.items():
        variant_to_refs[tuple(var)].append(header)
    amp_to_var[amp] = variant_to_refs

# add 




# check that removing variant doesn't make it indistinguable
#Counter(len(val) for val in ref_to_var.values()) 

# # min nb of unique variant per ref
# unique_amp = set()
# for amp,v_to_ref in amp_to_var.items():
#     for v,refs in v_to_ref.items():
#         if len(refs)==1:
#             unique_amp|={refs[0]}
# # check what's left
# ref_left = set(ref_to_var.keys())-unique_amp

# output variant combination
with open(OUT,"w") as handle:
    for amp,v_to_ref in amp_to_var.items():
        for v,refs in v_to_ref.items():
            if v==():
                continue
            variant = ";".join(["|".join(map(str,el)) for el in v])
            handle.write("%s\t%s\t%s\n"%(amp.replace(".msa",""),variant,"|delim|".join(refs)))
