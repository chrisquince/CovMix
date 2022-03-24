import os,glob
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp


for R1 in glob.glob("*L001_R1_001.fastq.gz"):
	sample = R1.split("_L001_R1_001.fastq.gz")[0]
	R2 = R1.replace("L001_R1_001.fastq.gz","L001_R2_001.fastq.gz")
	os.system("mkdir -p %s"%sample)
	os.system("mv %s* %s/"%(sample,sample))

for folder in glob.glob("*/"):
	os.system("mv %s %s"%(folder,folder[:-2]))




# first samples manipulation

barcode_to_sample = {line.rstrip().split("\t")[5]:line.rstrip().split("\t")[2] for line in open("metadata/barcode_to_sample.csv")}

fold_to_sample = {fold:barcode_to_sample[fold.replace("barcode0","").replace("barcode","")] for fold in  glob.glob("barcode*")}

for folder in glob.glob("barcode*"):
	sample = fold_to_sample[folder]
	os.system("mkdir -p %s/tmp"%folder)
	os.system("mv %s/*.fastq.gz %s/tmp/"%(folder,folder))
	os.system("cat %s/tmp/*.fastq.gz > %s/%s.fastq.gz &"%(folder,folder,sample)) 

for folder in glob.glob("barcode*"):
	sample = fold_to_sample[folder]
	os.system("mv %s %s"%(folder,sample))

header,seq = next(sfp(open("/mnt/gpfs/seb/Project/CovMix/primer_schemes/Artic_V4/Artic_V4_reference.fasta")))

# deal with primers
primers = []
with open("metadata/primers_suba_v46.csv") as handle:
	header = next(handle)
	for line in handle:
		name,pool,start,end,*_ = line.rstrip().split("\t")
		if name.split("_")[-1]=="F":
			side = "LEFT"
			orientation = "+"
		if name.split("_")[-1]=="R":
			side = "RIGHT"
			orientation = "-"
		pool = (pool=="e")*"2"+(pool=="o")*"1"
		# convert to 0 based
		new_start = int(start)-1
		new_end = int(end)-1
		new_name = "SARS-CoV-2_%s_%s"%(name.split("_")[2],side)
		primers.append(["MN908947.3",str(new_start),str(new_end),new_name,pool,orientation])
primers = sorted(primers,key=lambda x:int(x[3].split("_")[1]))

with open("primers_sub_4_6.bed","w") as handle:
	handle.writelines("%s\n"%"\t".join(line) for line in primers)


name_to_region2 = {}

suba_v46


