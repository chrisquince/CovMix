#!/usr/bin/env python3
import glob
import argparse
import numpy as np
from os.path import basename
import matplotlib.pyplot as plt
from operator import itemgetter
import matplotlib.backends.backend_pdf
from collections import defaultdict,Counter
from matplotlib.colors import ListedColormap
from Bio.SeqIO.FastaIO import SimpleFastaParser as sfp
# plot counts predicted vs observed and write that info in a file
# let's create a file with pos,base,percent


def union(intervals):
    sorted_intervals = sorted([tuple(sorted(i)) for i in intervals], key=itemgetter(0))
    if len(sorted_intervals)<=1:  # no intervals to merge
        return
    # low and high represent the bounds of the current run of merges
    low, high = sorted_intervals[0]
    for iv in sorted_intervals[1:]:
        if iv[0] <= high:  # new interval overlaps current run
            high = max(high, iv[1])  # merge with the current run
        else:  # current run is over
            yield low, high  # yield accumulated interval
            low, high = iv  # start new run
    yield low, high  # end the final run




# try another approach with completely unfiltered output
def unfiltered_proportions(ROOT,sorted_refs,sorted_samples):
    files = glob.glob("%s/samples/*/EM_runs/*_total_pi_est.csv"%ROOT)
    sample_to_file = {basename(file).replace("_total_pi_est.csv",""):file for file in files}

    # build a matrix
    mat = np.zeros((len(files),len(sorted_refs)))
    for index,sample in enumerate(sorted_samples):
        with open(sample_to_file[sample]) as handle:
            _=next(handle)
            for line in handle:
                _,ref,percent,std = line.rstrip().split(",")
                mat[index,sorted_refs.index(ref)]=float(percent)
    return mat,sorted_refs,sorted_samples

def get_coord(BED_FILE,amp_list):
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
    amp_to_coord = {key.split("_")[1]:[fwd[1],rev[0]] for key,(fwd,rev) in primer_to_coord.items()}
    ranges = [amp_to_coord[amp] for amp in amp_list]
    return list(union(ranges))

def degeneracy(mat,variants,sorted_sample):
    degen_summary = {(line.rstrip().split("\t")[0],line.rstrip().split("\t")[1]):line.rstrip().split("\t")[2].split("|----|") for index,line in enumerate(open("%s/results/degen_summary.tsv"%ROOT)) if index>0} 
    degen_summary = {(sample):[el for el in val if el!=rep] for (sample,rep),val in degen_summary.items()}
    variants_index = {var:index for index,var in enumerate(variants)}
    for index,sample in enumerate(sorted_sample):
        for variant in degen_summary[(sample)]:
            mat[index,variants_index[variant]]=0
    return mat


def load_matrix(file,sep="\t",sample_order=None,strain_order=None) :
    with open(file) as handle :
        header = next(handle).rstrip().split(sep)[1:]
        strains = []
        matrix = []
        for line in handle : 
            splitlines = line.rstrip().split(sep)
            strains.append(splitlines[0])
            matrix.append(list(map(float,splitlines[1:])))
    matrix = np.array(matrix)
    if sample_order :
        reorder_samples = [header.index(sample) for sample in sample_order]
        reorder_strain = [strains.index(strain) for strain in strain_order]
        return matrix[:,reorder_samples][reorder_strain,:]
    else : 
        return matrix,header,strains


def get_pos_base_percent(ROOT,prprt):
    # get snv 
    posbase_percent = defaultdict(lambda:np.zeros((len(sorted_sample))))
    for line in open("%s/refs/db_snv.tsv"%ROOT):
        variant,pos,base = line.rstrip().split("\t")
        base = base.upper()
        if (len(base)>1)|(base not in {"A","T","G","C"}):
            continue
        variant = variant.split(" ")[0]
        if variant not in variants:
            continue
        if sum(prprt[variants.index(variant),:])==0:
            continue
        pos = int(pos)
        posbase_percent[(pos,base)]+=prprt[variants.index(variant),:]
    return posbase_percent

def pos_base_map(ROOT):
    # get snv 
    snv_variants = defaultdict(list)
    for line in open("%s/refs/db_snv.tsv"%ROOT):
        variant,pos,base = line.rstrip().split("\t")
        base = base.upper()
        if (len(base)>1)|(base not in {"A","T","G","C"}):
            continue
        variant = variant.split(" ")[0]
        pos = int(pos)
        snv_variants[(pos,base)].append(variant)
    snv_variants = {key:"|sep|".join(val) for key,val in snv_variants.items()}
    return snv_variants


def get_observed_snv(ROOT,sorted_sample,varscan=True):
    if varscan:
        sample_to_pos = get_varscan_cnts(ROOT)
    cnts_matrix,sample_bases,ref = load_matrix("%s/snv/count.csv"%ROOT,sep=",")
    sample_to_indexes = defaultdict(list)
    for index,sample in enumerate(sample_bases[1:]):
        sample_to_indexes[sample[:-2]].append(index)
    sample_to_indexes = {key:np.array(val) for key,val in sample_to_indexes.items()}

    sorted_bases = ["A","C","G","T"]
    posbase_cnts_observed = defaultdict(lambda:np.zeros(len(sorted_sample)).astype(int))
    pos_cnts = defaultdict(lambda:np.zeros(len(sorted_sample)).astype(int))

    for index,line in enumerate(cnts_matrix):
        pos = int(line[0])-1
        line = line[1:]
        if sum(line)==0:
            continue
        for sample,indexes in sample_to_indexes.items():
            if varscan:
                if pos not in sample_to_pos[sample]:
                    continue
            tot_nuc = sum(line[indexes])
            pos_cnts[pos][sorted_sample.index(sample)] = int(tot_nuc)
            if tot_nuc!=0:
                for index,cnts in enumerate(line[indexes]):
                    base = sorted_bases[index]
                    if refseq[pos]!=base:
                        posbase_cnts_observed[(pos,base)][sorted_sample.index(sample)] = int(cnts)
    return posbase_cnts_observed,pos_cnts

def filter_combi(sample_combi,ranges):
    sorted_combi_init = sorted(sample_combi,key=lambda x:x[0])
    sorted_combi = []
    index_ranges = 0
    for pos,nuc in sorted_combi_init:
        if pos < ranges[index_ranges][0]:
            continue
        if pos > ranges[index_ranges][1]:
            if index_ranges<(len(ranges)+1):
                index_ranges+=1
            continue
        sorted_combi.append((pos,nuc))
    return sorted_combi

def get_varscan_cnts(ROOT):
    sample_to_pos = {}
    files = glob.glob("%s/samples/*/varscan/*_snv.tsv"%ROOT)
    for file in files:
        sample = basename(file).replace("_snv.tsv","")
        pos = {int(line.rstrip().split("\t")[1])-1 for index,line in enumerate(open(file)) if index>0}
        sample_to_pos[sample]=pos
    return sample_to_pos





def sample_wise_data(sample_index,posbase_cnts,posbase_cnts_observed,posbase_filt, sorted_sample, snv_variants, pos_cnts,ranges):
    sorted_cases = ["observed and predicted (%s)", "not predicted (%s)", "not observed (%s)", "observed only once (%s)", "removed by threshold (%s)"]
    sample_pred = {key for key,val in posbase_cnts.items() if val[sample_index]!=0}
    sample_obs = {key for key,val in posbase_cnts_observed.items() if val[sample_index]!=0}
    sample_combi = sample_pred|sample_obs
    sample_common = sample_pred&sample_obs
    color_label = {}
    sorted_combi = filter_combi(sample_combi,ranges)
    for combi in sorted_combi:
        if combi in sample_common:
            color_label[combi] = "observed and predicted (%s)" 
        if combi in sample_obs-sample_pred:
            color_label[combi] = "not predicted (%s)"
        if combi in sample_pred-sample_obs: 
            color_label[combi] = "not observed (%s)"
        if sum(np.array(posbase_cnts_observed[combi])!=0)<=1:
            color_label[combi] = "observed only once (%s)"
        if (combi in sample_common)|(combi in sample_pred-sample_obs):
            if combi in posbase_filt:
                if posbase_filt[combi][sample_index]==0:
                    color_label[combi] = "removed by threshold (%s)"
            if combi not in posbase_filt:
                color_label[combi] = "removed by threshold (%s)"

    label_to_val = {case:index for index,case in enumerate(sorted_cases)}
    color = [label_to_val[color_label[combi]] for combi in sorted_combi]
    nb = Counter(color_label.values())
    sample_cases = [case%nb[case] for case in sorted_cases]

    # get x,y
    x_predict = [posbase_cnts[comb][sample_index] if comb in posbase_cnts else 0 for comb in sorted_combi]
    y_obsv = [posbase_cnts_observed[comb][sample_index] if comb in posbase_cnts_observed else 0 for comb in sorted_combi]

    # hack: not all colors are shown if I don't do that;
    x_predict += [1e-10,1e-10,1e-10,1e-10,1e-10]
    y_obsv +=[1e-10,1e-10,1e-10,1e-10,1e-10]
    color +=list(range(len(sorted_cases)))

    # create data 
    data = []
    for combi in sorted_combi:
        pos,base = combi
        if (combi in snv_variants)&(combi in posbase_filt):
            origin = snv_variants[combi]
            count_pred = posbase_cnts[combi][sample_index]
        else:
            origin = "observed"
            count_pred = 0
        if combi in posbase_cnts_observed:
            count_obs = posbase_cnts_observed[combi][sample_index]
        else:
            count_obs = 0
        data.append([sorted_sample[sample_index],origin,pos,base,count_obs,count_pred,pos_cnts[combi[0]][sample_index]])
    return sample_cases, color, x_predict, y_obsv, data


parser = argparse.ArgumentParser(description='generate plot of predicted vs observed snv counts as well as flat file with the same info. Output in the snv folder of covmix')
parser.add_argument("folder", help='covmix folder path')
args = parser.parse_args()


ROOT = args.folder

# get reference file used for snv calling
ref_file = glob.glob("%s/refs/*_reference.fasta"%ROOT)[0]
name,refseq = next(sfp(open(ref_file)))

# get filtered predicted variant proportions
prprt_filt,variants,sorted_sample = load_matrix('%s/results/proportions.tsv'%ROOT)
prprt_filt = degeneracy(prprt_filt,variants,sorted_sample)
prprt_filt = prprt_filt.T

# get unfiltered variant proportions
prprt,variants,sorted_sample = unfiltered_proportions(ROOT,variants,sorted_sample)
prprt = degeneracy(prprt,variants,sorted_sample)
prprt = prprt.T


# get primer range
PRIMER_FILE = glob.glob("%s/refs/*_primer.bed"%ROOT)[0]
RANGE = [min([int(line.rstrip().split("\t")[1]) for line in open(PRIMER_FILE) if "LEFT" in line if "alt" not in line]),max([int(line.rstrip().split("\t")[1]) for line in open(PRIMER_FILE) if "LEFT" not in line if "alt" not in line])]

# for each sample assess region used for calling variants
sample_to_range = {}
for file in glob.glob("%s/samples/*/EM_runs/*_total_Filt_amp_pi_est.csv"%ROOT):
    sample = basename(file).replace("_total_Filt_amp_pi_est.csv","")
    amp_list = next(open(file)).rstrip().split(",")[2:]
    sample_to_range[sample]=get_coord(PRIMER_FILE,amp_list)

# get db snv mapping
snv_variants = pos_base_map(ROOT)


# get snv proportions
posbase_filt = get_pos_base_percent(ROOT,prprt_filt)
posbase_percent = get_pos_base_percent(ROOT,prprt)

# get observed counts
posbase_cnts_observed,pos_cnts = get_observed_snv(ROOT,sorted_sample)

issue = {combi for combi,values in  posbase_cnts_observed.items() if sum(pos_cnts[combi[0]]-values)>0}

# translate snv proportion in to countgs
posbase_cnts = {(pos,base):(freq*pos_cnts[pos]).astype(int) for (pos,base),freq in posbase_percent.items()}



# plot output
results = []
pdf = matplotlib.backends.backend_pdf.PdfPages("%s/snv/snv_plot.pdf"%ROOT)
for index_page in range(0,len(sorted_sample),3):
    fig = plt.figure(figsize=(10, 10))
    for nb_fig in range(3):
        plot_num = 311+nb_fig
        # plot corresponding dmags
        index_tot = index_page+nb_fig
        if (index_tot+1)>len(sorted_sample):
            continue
        plt.subplot(plot_num)

        # get data
        sample_cases, color, x_predict, y_obsv, data = sample_wise_data(index_tot,posbase_cnts,posbase_cnts_observed,posbase_filt,sorted_sample,snv_variants,pos_cnts,sample_to_range[sorted_sample[index_tot]])

        # add to result for futur tabular output
        results+=data

        # plot 1rst bisectrice
        _min,_max = int(min(x_predict+y_obsv)),int(max(y_obsv+x_predict))
        baseline = list(range(_min,_max,max(1,int((_max-_min)/100))))

        # sctter plot with pseudo count for log scale
        log_pseudo_count = 9e-1
        plt.plot(baseline,baseline,'b',alpha=0.5,linewidth=1.5)
        scatter = plt.scatter(np.array(x_predict)+log_pseudo_count,np.array(y_obsv)+log_pseudo_count,marker='o',alpha=0.5,c=color,cmap=ListedColormap(["g","m","c","k","r","y"]),s=20)
        classes= sample_cases
        scatter.legend_elements(prop='colors', num=len(sample_cases)) 
        # add legend outside of axes
        plt.legend(handles=scatter.legend_elements()[0], labels=classes,bbox_to_anchor=(1.04,0.5), loc="center left", borderaxespad=0)
        plt.xlabel("predicted snv counts")
        plt.ylabel("observed snv counts")
        plt.title("observed against predicted\n snv counts for sample %s"%(sorted_sample[index_tot]))
        plt.yscale('log')
        plt.xscale('log')
        plt.ylim(log_pseudo_count*0.9,1.5*_max)
        plt.xlim(log_pseudo_count*0.9,1.5*_max)
    fig.tight_layout()
    pdf.savefig(fig,bbox_inches="tight")
pdf.close()


# output a data table with the same info
with open("%s/snv/snv_table.tsv"%ROOT,"w") as handle:
    handle.write("%s\n"%"\t".join(["sample","origin","position","base","counts_observed","counts_predicted","total_counts"]))
    handle.writelines("%s\n"%"\t".join(map(str,line)) for line in results)




