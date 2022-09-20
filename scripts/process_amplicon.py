#!/usr/bin/env python3
from Bio.SeqIO.QualityIO import FastqGeneralIterator as fgi
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from common import *
import argparse
import gzip


def vsearch_mergepair(ROOT):
    R1 = f'{ROOT}_R1.fastq.gz'
    R2 = f'{ROOT}_R2.fastq.gz'
    output = f'{ROOT}.merged'
    cmd = f"vsearch --threads 1 --fastq_mergepairs {R1} --reverse {R2} --fastq_minovlen 5 --fastaout {output} || touch {output}"
    shell_cmd(cmd,"merged")

def deal_with_low_overlap(ROOT):
    R1_handle = fgi(gzip.open(f"{ROOT}_no_ovlp_R1.fastq.gz", "rt"))
    R2_handle = fgi(gzip.open(f"{ROOT}_no_ovlp_R2.fastq.gz", "rt"))
    with open(f'{ROOT}.todereplicate',"w") as handle:
        for h_R1,seq_R1,qual_R1 in R1_handle:
            h_R2,seq_R2,qual_R2 = next(R2_handle)
            assert(h_R1==h_R2),"issue with non overlaping reads order"
            header,overlap = h_R1.split("_ovlp=")
            overlap = int(overlap)
            if (len(seq_R1)>overlap)&(len(seq_R2)>overlap):
                if overlap==0:
                    handle.write(">%s\n%s\n"%(header,seq_R1+seq_R2))
                else:
                    seq = seq_R1[:-overlap]
                    for i in range(0,overlap):
                        nuc1 = seq_R1[-overlap+i]
                        qual1 = ord(qual_R1[-overlap+i])-33
                        nuc2 = seq_R2[i]
                        qual2 = ord(qual_R2[i])-33
                        nuc = max([(nuc1,qual1),(nuc2,qual2)],key=lambda x:x[1])
                        seq += nuc[0]
                    seq += seq_R2[overlap:]
                    handle.write(">%s\n%s\n"%(header,seq))
            else:
                seq = max([seq_R1,seq_R2],key=len)
                handle.write(">%s\n%s\n"%(header,seq))
        handle.write(open(f'{ROOT}.merged').read())
        shell_cmd(f"rm {ROOT}.merged","remove previous temporary file (.merged)")

def fastq_to_fasta(ROOT):
    cmd = f"seqtk seq -a {ROOT}.fastq.gz > {ROOT}.todereplicate"
    shell_cmd(cmd,"extract fasta from fastq")

def dereplication(ROOT,MIN_DEREPLICATED_CNT):
    cmd = f"vsearch --threads 1 --derep_fulllength {ROOT}.todereplicate --output {ROOT}.fasta -sizeout -minuniquesize {MIN_DEREPLICATED_CNT}"
    shell_cmd(cmd,"use vsearch for dereplication")
    
    # rm previous temp file
    cmd = f"rm {ROOT}.todereplicate"
    shell_cmd(cmd,"remove previous temporary file (.todereplicate)")


# ------ obsolete, this is built in in vsearch ------ 
def remove_singletons(ROOT,MIN_DEREPLICATED_CNT):
    with open(f'{ROOT}.tofilter') as handle, open(f'{ROOT}.fasta',"w") as handle_w:
        for header,seq in sfp(handle):
            nb = int(header.split()[0].split(";size=")[1])
            if nb<MIN_DEREPLICATED_CNT:
                continue
            handle_w.write(">%s\n%s\n"%(header,seq))
    cmd = f"rm {ROOT}.tofilter","remove previous temporary file (.tofilter)"

# ------ obsolete, we don't want to subsamples anymore from having derep ------ 
def subsample_reads(ROOT,MAX_READS):
    cmd = f"seqtk sample -s42 {ROOT}.fasta {MAX_READS} > {ROOT}_F{MAX_READS}.fasta"
    shell_cmd(cmd,"subsample files")
	# rm previous temp file
    cmd = f"rm {ROOT}.fasta"
    shell_cmd(cmd,"remove previous temporary file (.fasta)")

def vsearch_all_merged(ROOT,DB_AMP,MIN_PID):
    # remove gap and add dummy seq
    INPUT = f'{ROOT}.fasta'
    TMP = f'{ROOT}.tmp'
    with open(TMP,"w") as handle:
        handle.writelines(f'>{header}\n{seq.replace("-","")}\n' for header,seq in sfp(open(DB_AMP)))
        handle.write('>need_at_least_one_entry_so_vsearch_doesnt_throw_an_error\nATGC\n')

    # vsearch
    cmd = f'vsearch --threads 1 --usearch_global {INPUT} --db {TMP} --samout {ROOT}.sam --samheader --id {MIN_PID} --maxaccepts 1000000 --userout {ROOT}.m6 --userfields query+target+evalue+id+pctgaps+pairs+gaps+qlo+qhi+tlo+thi+pv+ql+tl+qs+ts+alnlen+opens+exts+raw+bits+aln+caln+qrow+trow+mism+ids+qcov+tcov'
    shell_cmd(cmd,"run vsearch global alignment")

    # remove temp files
    cmd = f"rm {TMP}"
    shell_cmd(cmd,"remove previous temporary file (DB_AMP.tmp)")


def process_amplicon(ROOT, DATATYPE, MIN_DEREPLICATED_CNT, DB_AMP, MIN_PID):
    if DATATYPE=="p-reads":
        vsearch_mergepair(ROOT)
        deal_with_low_overlap(ROOT)
    else:
        fastq_to_fasta(ROOT)
    dereplication(ROOT, MIN_DEREPLICATED_CNT)
    vsearch_all_merged(ROOT, DB_AMP, MIN_PID)


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("prefix",help="prefix to R1.fastq.gz/R2.fastq.gz amplicon file")
    parser.add_argument("datatype",help="is this paired reads or nanopore?")
    parser.add_argument("min_dereplicated_cnt",help="threshold for minimum cnt an amplicon has been seen")
    parser.add_argument("amp_db",help="amplicon specifique database")
    parser.add_argument("min_pid",help="minimum identity for global alignment")
    args = parser.parse_args()

    process_amplicon(args.prefix, args.datatype, args.min_dereplicated_cnt, args.amp_db, args.min_pid)
    

