import sys
import os
import argparse
import glob
from ete3 import NCBITaxa


import clusteread
import taxtree

# remove all files inside a folder
def rmf(paths):
    for path in paths:
        for file in glob.glob(os.path.join(path, "/*")):
            os.remove(file)

def clustering(path):
    pass

# reqs: id/
def identity(path):
    pass

# reqs: hmmsearch/
def similarity(path):
    pass

def find_highest_id(ncbi, fasta):
    id2taxid = taxtree.read_cluster_taxids(fasta)
    tax_tree = ncbi.get_topology(id2taxid.values(), intermediate_nodes=False)
    highest_tax = taxtree.get_highest_taxonomic_id(id2taxid.values(), tax_tree)
    for id in id2taxid.keys():
        if id2taxid[id] == highest_tax:
            return id
    
    assert False, "Highest taxid not found..."

def set_ref(path, ref):
    data = {}
    with open(path, "r") as f:
        db_fasta = ("\n" + f.read()).split("\n>")[1:]
        for protein in db_fasta:
            protein_id, sequence = clusteread.parse_fasta(protein)
            clean_id = protein_id.split("|")[1]
            data[clean_id] = protein
    
    with open(path, "w") as f:
        f.write(">" + data[ref] + "\n")
        for id in data.keys():
            if id == ref:
                continue
            f.write(">" + data[id] + "\n")
    
    if path.contains("clean"):
        ext = ".clean.fasta"
    else:
        ext = ".fasta"
    os.rename(path, os.path.join("".join(path.split("/")[:-1]), ref, ext))

def taxonomy(path):
    rmf([path + "/phmmer", path + "/refs"])
    fasta_files = glob.glob(os.path.join(path, "/fasta/*"))
    clean_files = glob.glob(os.path.join(path, "/clean/*"))
    files = zip(fasta_files, clean_files)
    ncbi = NCBITaxa()
    for file in files:
        reference = find_highest_id(ncbi, file[0])
        set_ref(path[0], reference)
        set_ref(path[1], reference)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    # parser.add_argument("db", type=str, help="Database in fasta-format")
    parser.add_argument("path", type=str, help="Path to the cluster group directory")
    # parser.add_argument("criterion", type=str, help="reference selection criterion: d(efault: first node in multi-step/mcl)/i(dentity)/h(mmsearch score)")
    parser.add_argument("--clustering", action="store_true", help="First node in multi-step/mcl")
    parser.add_argument("--identity", action="store_true", help="Highest mean pair-wise identity score")
    parser.add_argument("--similarity", action="store_true", help="Highest hmmsearch score (i.e. most similar sequence to the MSA of the cluster)")
    parser.add_argument("--taxonomy", action="store_true", help="Highest node in taxonomic tree")
    args = parser.parse_args()
    if args.clustering:
        clustering(args.path)
    elif args.identity:
        identity(args.path)
    elif args.similarity:
        similarity(args.path)
    elif args.taxonomy:
        taxonomy(args.path)
