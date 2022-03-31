#!/usr/bin/env python3
import argparse
import pandas as pd
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


parser = argparse.ArgumentParser()
parser.add_argument("metadata",help="GISAID metadata.tsv file, should be the uncompressed flat file")
parser.add_argument("seq",help="GISAID sequences.fasta file of all genomes, should be the uncompressed flat file")
parser.add_argument("outfold",help="output folder where representatives.fa and rep_metadata.tsv will be written")


args = parser.parse_args()  
metadata = args.metadata
output_folder = args.outfold
seq_file = args.seq


df = pd.read_csv(metadata,sep='\t')

# make header normals
df.columns =  df.columns.map(lambda x:x.lower().replace(" ","_").replace("-","_").replace("?",""))

# filter for human host and complete and high coverage
df = df.query("host == 'Human' and is_complete == True and is_high_coverage == True")

# for each pango_lineage find the representative with the least nb of N
selected = df.groupby('pango_lineage')["n_content"].nsmallest(1).reset_index().iloc[:,1]
rep_df = df.loc[selected,]

# output infos 
rep_df.to_csv("%s/rep_metadata.tsv"%output_folder, sep='\t', index=False)

ids_lineages = rep_df.set_index('virus_name')['pango_lineage'].to_dict()
rep_ids = set(ids_lineages.keys())

with open("%s/representatives.fa"%output_folder,"w") as handle:
    for header,seq in sfp(open(seq_file)):
        name = header.split("|")[0]
        if name in rep_ids:
            lineage = ids_lineages[name]
            handle.write(">%s|%s\n%s\n"%(name,lineage,seq))

