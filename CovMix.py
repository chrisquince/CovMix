#!/usr/bin/env python3
from os.path import abspath, realpath, dirname, basename, exists
from scripts.common import fill_default_values, cd 
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
parser.add_argument("--dag", "-d", help="file where you want the dag to be stored")
parser.add_argument('-s', nargs=argparse.REMAINDER,help="Pass additional argument directly to snakemake")
args = parser.parse_args()

# get config file
CONFIG_FILE = abspath(realpath(args.config))
config = yaml.full_load(open(CONFIG_FILE))

# get repos directory
REPOS_DIR = dirname(abspath(realpath(sys.argv[0])))

# execution directory
EXEC_DIR=abspath(realpath(config["execution_directory"]))
os.system("mkdir -p %s"%EXEC_DIR)

# ------- base parameters used to call snakemake -----------
base_params = ["snakemake", "--directory", EXEC_DIR, "--cores", str(args.cores), "--config", "REPOS_DIR=%s"%REPOS_DIR,"PRIMER=%s"%args.primer,"CONFIG_PATH=%s"%CONFIG_FILE,"EXEC_DIR=%s"%EXEC_DIR,"DATATYPE=%s"%args.datatype,"--configfile="+CONFIG_FILE, "--latency-wait", "120","-k","--use-conda"]

# ------- additional parameters -----------
if args.verbose:
    base_params.extend(["-p", "-r", "--verbose"]) 
if args.dryrun:
    base_params.extend(["--dryrun"])
if args.unlock:
    base_params.extend(["--unlock"])
if args.touch:
    base_params.extend(["-t"])    
if args.dag:
    base_params.extend(["--rulegraph"])
if args.s :
    base_params.extend(args.s)


# ------- call snakemake from  METAHOOD_DIR -----------
with cd(REPOS_DIR):
    def call_snake(extra_params=[]):
        call_snake.nb+=1
        if args.dag:
            p1=Popen(base_params + extra_params, stdout=PIPE, stderr=sys.stderr)
            p2=Popen(["dot","-Tpng"],stdin=p1.stdout, stdout=PIPE, stderr=sys.stderr)
            with open(args.dag.replace(".png",str(call_snake.nb)+".png"),"bw") as f :
                f.write(p2.communicate()[0])
        else :
            subprocess.check_call(base_params + extra_params, stdout=sys.stdout, stderr=sys.stderr)
    call_snake.nb=0
    # so to remove checkpoint the pipeline is separated in 2 bits: 
    # first trim and select amp to run
    call_snake(["--snakefile", "scripts/main_pipe.snake" ,"%s/selected_amp.tsv"%EXEC_DIR])
    # Then run regular EM, snv as well as varscan
    call_snake(["--snakefile", "scripts/main_pipe.snake" ])













