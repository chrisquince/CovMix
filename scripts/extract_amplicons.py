#!/usr/bin/env python3
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from collections import defaultdict,Counter
from multiprocessing import Pool
import numpy as np
import argparse
import glob
import time
import os
import re



parser = argparse.ArgumentParser(description='Extract amplicons using primer, if primer region is mutated, amplicon is missing so to mimic expected output')
parser.add_argument('primer_bed', help='bed file describing primer coordinates, beware this script assume some naming conventions, (_LEFT,_RIGHT)')
parser.add_argument('ref', help='fasta file with reference fasta used for primer definition')
parser.add_argument('rep', help='fasta file of genome databases')
parser.add_argument('amp', help='prefix for amplicon database sequences')
parser.add_argument('--next_of_kin', help="custom tsv file generated by tabulate_relatives mapping each leave of the tree to all (leave:distance), this file will be used to replace any missing amplicon by the nearest relative's one" )
parser.add_argument('-m', help='minimum number of clean amplicon for a ref to be kept, depend on database/primer scheme')
parser.add_argument('-s', help='output some stats on database: nb unique amplicons/max nb of shared amplicons ....')
args = parser.parse_args()

# REF = "/mnt/gpfs/seb/Database/covid/nCoV-2019.reference.fasta"
# BED_FILE = "/mnt/gpfs/seb/Database/covid/nCoV-2019.primer.bed"
# REP = "/mnt/gpfs/seb/Database/covid/GISAID/reps_and_ref.fa"
# TREE = "/mnt/gpfs/seb/Database/covid/GISAID/reps_and_ref.tree"
# TREE_REL = "/mnt/gpfs/seb/Database/covid/GISAID/reps_and_ref_tree_rel.tsv"

# --------- custom functions --------- 
# use primer to fish out amplicons
def get_ampseq(amp_to_primer,seq):
    ampseq = {}
    missing = {}
    degenerate = {}
    for amp,(p1,p2) in amp_to_primer.items():
        m = re.findall("%s(.*)%s"%(p1,p2),seq.upper().replace("-",""))
        if m:
            if set(m[0])-{"A","T","G","C"} == set():
                ampseq[amp]=m[0]
            else:
                degenerate[amp]=""
        else:
            missing[amp]=""
    return ampseq,missing,degenerate



# ------------ start of script ---------------

# pass args, maybe useless
BED_FILE = args.primer_bed
REF = args.ref
REP = args.rep
AMP_OUT = args.amp

# filter refs by min number of amplicons
MIN_AMP = 0
if args.m : 
    MIN_AMP = int(args.m)

# replace missing amplicons (degenerated base of primer not found) by nearest ref version
if args.next_of_kin:
    TREE_REL = args.next_of_kin
else:
    TREE_REL =""

# output some stats on relatedness between refferences.
if args.s:
    STAT_OUT = args.s
else:
    STAT_OUT = AMP_OUT


# -------- real start of script --------
# get refference data
header,seq_ref = next(sfp(open(REF)))
amp_to_coord = {line.rstrip().split("\t")[3]:list(map(int,line.rstrip().split("\t")[1:3])) for line in open(BED_FILE)}

# get primer sequences
amp_to_primer = defaultdict(lambda:[None,None])
for amp_primer,coord in amp_to_coord.items():
    amp = amp_primer.split("_RIGHT")[0].split("_LEFT")[0]
    amp_to_primer[amp]["RIGHT" in amp_primer]=seq_ref[coord[0]:coord[1]]


# get representatives amplicons sequences
header_to_ampseq = {}
header_to_miss = {}
header_to_degenerate = {}
for header,seq in sfp(open(REP)):
    header = header.split(" ")[0]
    header_to_ampseq[header],header_to_miss[header],header_to_degenerate[header] = get_ampseq(amp_to_primer,seq)

header_to_fail = {header:{amp:"" for amp in set(ampseq.keys())|set(header_to_degenerate[header].keys())} for header,ampseq in header_to_miss.items()}

header_to_miss = {header:ampseq for header,ampseq in header_to_miss.items() if ampseq!={}}
header_to_degenerate = {header:ampseq for header,ampseq in header_to_degenerate.items() if ampseq!={}}


# fix degenerate|nonfound amplicon sequences
if TREE_REL:
    header_to_rel = {line.rstrip().split("\t")[0]:[el.split(":")[0] for el in line.rstrip().split("\t")[1:] if el.split(':')[0]!=line.rstrip().split("\t")[0]] for line in open(TREE_REL)}
    for header,ampseq in header_to_fail.items():
        for amp in ampseq:
            seq = next(header_to_ampseq[rel][amp] for rel in header_to_rel[header] if amp in header_to_ampseq[rel])
            header_to_ampseq[header][amp]=seq


# # --------- generate some stats --------- 
# # get the number of missing amp per amp
# degenerate = sorted(Counter([amp for ampseq in header_to_degenerate.values() for amp in ampseq]).items(),key=lambda x:x[1])
# missing = sorted(Counter([amp for ampseq in header_to_miss.values() for amp in ampseq]).items(),key=lambda x:x[1])
# # let's do something bruteforce : test removal of all possible amp. And find the best set of removed amplicon so that we keep the maximum number of refs.

# # create a handy matrix, ref against ref against amp presence so that we can do some fast matrix sum, to know what happens in terms of ambiguity when removing an amplicon
# amplicons = {amp for ampseq in header_to_ampseq.values() for amp in ampseq}
# ambiguity_mat = np.zeros((len(header_to_ampseq),len(header_to_ampseq),len(amplicons)),dtype=int)
# header_to_index = {header:ind for ind,header in enumerate(sorted(header_to_ampseq.keys()))}
# amp_to_index = {header:ind for ind,header in enumerate(sorted(amplicons))}

# for header,ampseq in header_to_ampseq.items():
#     for amp in ampseq:
#         for h in new_name_to_ref[ampheader_name[(amp,header)]]:
#             if header!=h:
#                 ambiguity_mat[header_to_index[header],header_to_index[h],amp_to_index[amp]]=1
# # create a ref x amp matrix showing nb of present amplicons
# completion_mat = np.zeros((len(header_to_ampseq),len(amplicons)),dtype=int)
# for header,ampseq in header_to_ampseq.items():
#     for amp in ampseq:
#         completion_mat[header_to_index[header],amp_to_index[amp]]=1


# # not optimal but I don't want to test all combinations, so I select the best option at n amplicons removed and use the best result as seed for n+1
# def get_stats(args):
#     amp,candidate_amps = args
#     indexes = np.array([amp_to_index[a] for a in candidate_amps if a!=amp])
#     L = len(indexes)
#     complete = (completion_mat[:,indexes].sum(1)==L)
#     unambiguous = ambiguity_mat[:,:,indexes].sum(2).max(1)!=L
#     return [amp,sum(complete),sum(unambiguous),sum(complete*unambiguous)]

# def recursive_removal(set_removed,nb_to_stat):
#     start = time.time()
#     candidate_amps = set(amp_to_index.keys())-set_removed
#     pool = Pool(min(THREADS,len(candidate_amps)))
#     args = [(amp,candidate_amps) for amp in candidate_amps]
#     res = pool.map(get_stats,args) 
#     best = max(res,key=lambda x:x[3])
#     nb_to_stat[len(set_removed)+1] = best
#     set_removed|={best[0]}
#     print("{0:4g}".format((time.time()-start)/60.),best)
#     if len(set_removed)==len(amp_to_index):
#         return
#     else:
#         recursive_removal(set_removed,nb_to_stat)


# # this delete 1 of the amplicon which impact less ambiguity and completion recursively 
# nb_to_stat = {}
# %time recursive_removal(set(),nb_to_stat)









# header_to_shared[header] = max(Counter([h for amp in ampseq for h in new_name_to_ref[ampheader_name[(amp,header)]] if h!=header]).items(),key=lambda x:x[1])


# def unambiguous_refs(to_ignore):
#     ampheader_name[(amp,header)

# max(Counter([h for amp in ampseq for h in new_name_to_ref[ampheader_name[(amp,header)]] if amp not in to_ignore if h!=header]).items(),key=lambda x:x[1])

















# filter amplicon sequence and only keep the one with at least MIN_AMP seq 
header_to_ampseq = {header:ampseq for header,ampseq in header_to_ampseq.items() if len(ampseq)>=MIN_AMP}

# get reference seq
ref_ampseq,_,_ = get_ampseq(amp_to_primer,seq_ref)

# output amp seqs
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










