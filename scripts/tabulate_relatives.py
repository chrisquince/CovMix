#!/usr/bin/env python3
from multiprocessing import Pool
from ete3 import Tree
import numpy as np
import argparse


parser = argparse.ArgumentParser(description='go through a phylogenetic tree and use ete3 to output distances between leaves')
parser.add_argument('tree', help='path to tree')
parser.add_argument('out', help='path to output')
parser.add_argument('-t',default=1, help='nb threads, can be quite slow without a few threads')
args = parser.parse_args()


# go through the tree and find the closest relatives
def get_ordered_relatives(args):
    t,leave_ref = args
    leaves = np.array(sorted(t.get_leaves(),key=lambda x:x.name))
    dists = [t.get_distance(leave_ref, l, topology_only=False) for l in leaves]
    new_order = np.argsort(dists)
    return [leave_ref.name,list(np.array(["%s:%s"%(l.name,"{:.4g}".format(dists[index])) for index,l in enumerate(leaves)])[new_order])]


# ------------ start of script ---------------

# pass args, maybe useless
TREE = args.tree
OUT = args.out
THREADS = int(args.t)

# go iterate trhough the tree and just get all the distances between nodes
# does all distances twices.... nobrainer this way
t = Tree(TREE)
pool = Pool(THREADS)
args = [(t,leave) for leave in t.get_leaves()]
ref_to_relatives = {refname:nearest_rel for refname,nearest_rel in pool.map(get_ordered_relatives,args)}

# output results
with open(OUT,"w") as handle:
	handle.writelines("%s\t%s\n"%(refname,"\t".join(rel)) for refname,rel in ref_to_relatives.items())

