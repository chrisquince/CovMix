from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
from Bio.Seq import Seq

name,refseq = next(sfp(open("Nimagen_V3_reference.fasta")))

# generate Nimagen_V3_primer.bed
with open("Nimagen_V3_primer.bed","w") as handle:
    for index,line in enumerate(open("SARS-CoV-2_Nimagen-v7_kit-V3.bed")):
        ref,start,end,seq1,seq2,pool = line.rstrip().split("\t")
        # coordinate are written in zero based, but also end coordinate is inclusive and can't be used in a python range : 
        # correct used of these specific coordinate is range(start-1,end) 
        # translate to 0 based,
        start = int(start)-1
        end = int(end)
        handle.write("%s\t%s\t%s\tnCoV-2019_%s_LEFT\t%s\t+\n"%(ref,start-len(seq1),start,index+1,pool.split("_")[0]))
        handle.write("%s\t%s\t%s\tnCoV-2019_%s_RIGHT\t%s\t-\n"%(ref,end,end+len(seq2),index+1,pool.split("_")[0]))


# amplicon this time
with open("/mnt/gpfs/seb/Project/CovMix/primer_schemes/Nimagen_V3/Nimagen_V3_amplicons.bed","w") as handle:
    for index,line in enumerate(open("/mnt/gpfs/seb/Project/CovMix/primer_schemes/Nimagen_V3/SARS-CoV-2_Nimagen-v7_kit-V3.bed")):
        ref,start,end,seq1,seq2,pool = line.rstrip().split("\t")
        # coordinate are written in zero based, but also end coordinate is inclusive and can't be used in a python range : 
        # correct used of these specific coordinate is range(start-1,end) 
        # translate to 0 based,
        start = int(start)-1
        end = int(end)
        handle.write("%s\t%s\t%s\tnCoV-2019_%s\n"%(ref,start,end,index+1))



# generate Nimagen_V2_primer.bed
with open("Nimagen_V2_primer.bed","w") as handle:
    for index,line in enumerate(open("primers_bedpe_MN908947.3V2.bed")):
        ref,start1,end1,ref2,start2,end2,name = line.rstrip().split("\t")
        # coordinate are fine here
        handle.write("%s\t%s\t%s\t%s_LEFT\t/\t+\n"%(ref,start1,end1,name))
        handle.write("%s\t%s\t%s\t%s_RIGHT\t/\t-\n"%(ref,start2,end2,name))


# generate Artic_V4 schemes
# dl files from https://github.com/artic-network/artic-ncov2019/tree/master/primer_schemes/nCoV-2019/V4
# just change names from SARS-CoV-2.reference.fasta to Artic_V4_reference.fasta