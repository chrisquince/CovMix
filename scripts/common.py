from __future__ import print_function
from __future__ import division

try:
    from future_builtins import zip
except:
    pass

from collections import defaultdict
from subprocess import Popen, PIPE
from os.path import basename, dirname, isfile
import os.path
import sys
import os
import re


default_values = {"trimming_strictness":1,"Proportion_Threshold":0.025, "data_regex":["*"],"threads":20,"min_reads_per_amp":40,"max_reads_per_amp":5000,"amplicon_to_run":"","genome_to_run":""}
default_values_v2 = lambda path,primer:{
"database":
    {
    "tree":"%s/refs/tree/db_tree.nwk"%path,
    "tree_rel":"%s/refs/tree/db_tree_rel.tsv"%path,
    "amplicons":"%s/refs/%s_amps_seqs"%(path,primer)
    }
} 


# ---- neat regex matching of files --------
def extended_glob(pattern):
    process = Popen(['bash -O extglob -c " ls -d '+pattern+'/ "'], stdout=PIPE, stderr=PIPE,shell=True)
    List_path=[element[:-1] for element in process.communicate()[0].decode("utf-8").split("\n") if element]
    return [path for path in List_path if basename(path)!="multiqc_data"]


# Taken from http://stackoverflow.com/questions/36831998/how-to-fill-default-parameters-in-yaml-file-using-python
def setdefault_recursively(tgt, default = default_values):
    for k in default:
        if isinstance(default[k], dict): # if the current item is a dict,
            # expand it recursively
            setdefault_recursively(tgt.setdefault(k, {}), default[k])
        else:
            # ... otherwise simply set a default value if it's not set before
            tgt.setdefault(k, default[k])

def fill_default_values(config):
    REPOS_DIR = config.get("REPOS_DIR")
    EXEC_DIR = config.get("EXEC_DIR")
    PRIMER = config.get("PRIMER")
    default_values["scripts"] = os.path.join(REPOS_DIR, "scripts")
    setdefault_recursively(config)
    setdefault_recursively(config,default=default_values_v2(EXEC_DIR,PRIMER))


FASTA_EXTS = {".fastq.gz",".fastq",".fq.gz",".fq",".fa",".fasta",".fasta.gz",".fa.gz"}  # only extension valid. 
def gather_paths(path):
    for filename in os.listdir(path):
        name = os.path.basename(filename)
        for ext in FASTA_EXTS:
            if not name.endswith(ext):
                continue
            if "trimmed" in name : 
                continue
            if "Filtered" in name :
                continue
            yield os.path.join(path, filename)



def detect_reads(dir):
    return sorted(list(gather_paths(dir)))

def get_extension(file) :
    for ext in FASTA_EXTS:
        if file.endswith(ext):
            return ext 

def replace_extensions(sample,FILTER):
    if FILTER:
        sample = "%s/Filtered_%s"%(dirname(sample),basename(sample))
    ext = get_extension(sample)
    if ext in {".fastq.gz",".fastq",".fq.gz",".fq"} :
        return sample.replace(ext,"_trimmed%s"%ext)
    else :
        return sample

#copied from http://stackoverflow.com/questions/431684/how-do-i-cd-in-python/13197763#13197763
class cd:
    """Context manager for changing the current working directory"""
    def __init__(self, newPath):
        self.newPath = os.path.expanduser(newPath)

    def __enter__(self):
        self.savedPath = os.getcwd()
        os.chdir(self.newPath)

    def __exit__(self, etype, value, traceback):
        os.chdir(self.savedPath)

# assess samples :
def assess_samples(SAMPLES):
    Files_nb = []
    R = []
    Files_empty = set()
    for sample,files in SAMPLES.items():
        check = "%s/sample_ok"%dirname(files[0])
        if isfile(check):
            continue
        # check that there is exactly 2 samples
        if len(files)!=2:
            Files_nb.append(sample)
        # check that there is R1 and R2
        R1 = sum("_R1" in file for file in files)==0
        R2 = sum("_R2" in file for file in files)==0
        if R1|R2:
            R.append(sample)
        # check that sample is not empty
        for file in files:
            command = ["zcat","-f",file]
            p1 = Popen(command, stdout=PIPE, stderr=sys.stderr)
            p2 = Popen(["head"],stdin=p1.stdout, stdout=PIPE, stderr=sys.stderr)
            if len(p2.communicate()[0].split(b"\n"))<10:
                Files_empty.add(sample)
        if (not Files_empty)&(not R)&(not Files_nb):
            os.system("touch %s"%check)
    if Files_nb:
        sys.exit("Samples folder : "+" - ".join(Files_nb)+"  do not have exactly 2 reads files, you may have more or less than 2, files recognised as reads files are the following : .fastq, .fastq.gz, .fq, .fq.gz, .fa, .fa.gz, .fasta, .fasta.gz. You may also want to check the regular expression you used to select samples")
    if R:
        sys.exit("Samples folder : "+" - ".join(R)+" is missing _R1/_R2 tags in files names, please correct that")
    if Files_empty:
        sys.exit("Samples folder : "+" - ".join(Files_empty)+" do have R1 or R2 file with less than 10 lines, please remove these folder from analysis")
    if len(SAMPLES)==0:
        sys.exit("No sample folder has been detected, please ensure the path corresponding to data lead to a folder containing 1 folder per sample.")

# assess samples :
def assess_nanopore_samples(SAMPLES):
    Files_nb = []
    R = []
    Files_empty = set()
    for sample,files in SAMPLES.items():
        check = "%s/sample_ok"%dirname(files[0])
        if isfile(check):
            continue
        # check that there is exactly 1 fastq
        if len(files)!=1:
            Files_nb.append(sample)
        # check that sample is not empty
        for file in files:
            command = ["zcat","-f",file]
            p1 = Popen(command, stdout=PIPE, stderr=sys.stderr)
            p2 = Popen(["head"],stdin=p1.stdout, stdout=PIPE, stderr=sys.stderr)
            if len(p2.communicate()[0].split(b"\n"))<10:
                Files_empty.add(sample)
        if (not Files_empty)&(not Files_nb):
            os.system("touch %s"%check)
    if Files_nb:
        sys.exit("Samples folder: "+" - ".join(Files_nb)+"  do not have exactly 1 reads file, you may have more or less than 1, files recognised as reads files are the following : .fastq, .fastq.gz, .fq, .fq.gz, .fa, .fa.gz, .fasta, .fasta.gz. You may also want to check the regular expression you used to select sample folders")
    if Files_empty:
        sys.exit("Samples folder: "+" - ".join(Files_empty)+" is less than 10 lines long, maybe empty, please remove these folder from analysis")
    if len(SAMPLES)==0:
        sys.exit("No sample folder has been detected, please ensure the path corresponding to data lead to a folder containing 1 folder per sample.")