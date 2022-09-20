#!/usr/bin/env python3
from os.path import abspath, realpath, dirname, basename, exists
from scripts.common import *
from subprocess import PIPE,Popen
from psutil import virtual_memory
import subprocess
import argparse
import shutil
import yaml
import time
import sys
import os



parser = argparse.ArgumentParser(description="CovMix - pipeline estimating variant proportions from amplicons")
parser.add_argument('primer',
                    choices=['Artic_V3', "Artic_V4", "Artic_V4-6" , 'Nimagen_V2', 'Nimagen_V3', "Nimagen_V4"],
                    help='specify primer scheme, from 5 currently supported : Artic_V4, Artic_V3, Nimagen_V2, Nimagen_V3, Nimagen_V4, Artic_V4-6')
parser.add_argument("config", type=str, help="config_file.yaml to use")
parser.add_argument("--cores", "-c", type=int, default=1, help="Number of threads")
parser.add_argument('--datatype', choices=["p-reads","ont"],default="p-reads",help="Illumina short paired reads or nanopore long reads")
parser.add_argument("--verbose", "-v", action="store_true", help="Increase verbosity level")
parser.add_argument("--dryrun", "-n", action="store_true", help="Show tasks, do not execute them")
parser.add_argument("--unlock", "-u", action="store_true", help="Unlock the directory")
parser.add_argument("--touch", "-t", action="store_true", help="Touch all files, to reset timestamp and stop unwanted/uncalled snakemake reruns")
parser.add_argument("--dag", "-d", help="file where you want the dag to be stored",default="")
parser.add_argument("--mode", "-m", help="optionally run a unique part of the pipeline",choices=["all","map","amplicons","em","snv"],default="all")
parser.add_argument('-s', nargs=argparse.REMAINDER,help="Pass additional argument directly to snakemake")
args = parser.parse_args()

# get config file
CONFIG_FILE = abspath(realpath(args.config))
config = yaml.full_load(open(CONFIG_FILE))

# get repos directory
REPOS_DIR = dirname(abspath(realpath(sys.argv[0])))

# execution directory
EXEC_DIR = abspath(realpath(config["execution_directory"]))
os.system("mkdir -p %s"%EXEC_DIR)

# ------- base parameters used to call snakemake -----------
base_params = ["snakemake", "--directory", EXEC_DIR, "--cores", str(args.cores), "--configfile="+CONFIG_FILE, "--latency-wait", "120","-k","--use-conda"]
config_args = ["--config", "REPOS_DIR=%s"%REPOS_DIR,"PRIMER=%s"%args.primer,"CONFIG_PATH=%s"%CONFIG_FILE,"EXEC_DIR=%s"%EXEC_DIR,"DATATYPE=%s"%args.datatype]
# ------- additional parameters -----------
if args.verbose:
    base_params.extend(["-p", "-r", "--verbose"]) 
if args.dryrun:
    base_params.extend(["--dryrun"])
if args.unlock:
    base_params.extend(["--unlock"])
if args.touch:
    base_params.extend(["-t"])    
if args.s :
    base_params.extend(args.s)
    NOTHING_ELSE = 1
else:
    NOTHING_ELSE = 0

# -------- Samples --------
SAMPLE_DIR = config["data"]
try : 
    REGEX = config["data_regex"]
except:
    REGEX = ["*"]

SAMPLES = {basename(file):detect_reads(file) for regex in REGEX for file in extended_glob("%s/%s"%(SAMPLE_DIR,regex))}
SAMPLES = {sample:files for sample,files in SAMPLES.items() if files}

# check there is only 1 R1 and 1 R2 per folder
if args.datatype=="p-reads":
    assess_samples(SAMPLES)
    R1 = {sample:files["_R1" in files[1]] for sample,files in SAMPLES.items()}
    R2 = {sample:files["_R2" in files[1]] for sample,files in SAMPLES.items()}
else: 
    assess_nanopore_samples(SAMPLES)
    SAMPLES = {sample:files[0] for sample,files in SAMPLES.items()}

# ---------------------------------------------------------------
# ---------------------    Main pipeline    ---------------------
# ---------------------------------------------------------------
mode = args.mode
# cd allow to work in the correct dir 
with cd(REPOS_DIR):
    call_snake.nb=0 # thing for dag pic generation
    
    # start by parallel mapping all samples to ref
    if (mode=="all")|(mode=="map"):
        call_snake(base_params+["--snakefile", "scripts/map_to_ref.snake"]+config_args,"Mapping samples to database")
    
    # go back to sample wise processing
    if (mode=="all")|(mode=="amplicons"):
        for index,sample in enumerate(SAMPLES):
            call_snake(base_params+["--snakefile", "scripts/vsearch_amplicons.snake"]+config_args+["SAMPLE=%s"%sample],"Amplicon-wise global alignemnt for sample %s (%s/%s)"%(sample, index+1, len(SAMPLES)))
    
    # assess relevant amp/sample, over all samples
    if (mode=="all")|(mode=="em"):
        call_snake(base_params+["--snakefile", "scripts/EM_results.snake"]+[f"{EXEC_DIR}/selected_amp.tsv"]+config_args,"Assessing amplicons and samples to run")

        # recall the same snake but this time for launching all EM and such.
        call_snake(base_params+["--snakefile", "scripts/EM_results.snake"]+config_args,"Running EM algorithm, it may take some times")

    # do some additional snv analysis
    if (mode=="all")|(mode=="snv"):
        call_snake(base_params+["--snakefile", "scripts/snv.snake"]+config_args,"SNV analysis")

    print('-'*80)
    print("\n\n                  We are done here,\n\n                  ~ bye now\n\n")



