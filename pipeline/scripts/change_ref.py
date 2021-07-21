import sys
import os
import argparse
import glob
from ete3 import NCBITaxa
import pandas as pd

import clusteread
import taxtree

# remove all files inside a folder
def rmf(paths):
    for path in paths:
        files = glob.glob(os.path.join(path, "*"))
        for file in files:
            os.remove(file)

def clean_column_ids(df, col):
    return df[col].map(lambda x: x.split("|")[1])

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

    if "clean" in path:
        ext = ".clean.fasta"
    else:
        ext = ".fasta"
    os.rename(path, os.path.join("/".join(path.split("/")[:-1]), ref + ext))
    return data[ref]

def set_refs(path, ref, file_zip):
    set_ref(file_zip[0], ref)
    ref_fasta = set_ref(file_zip[1], ref)
    id, seq = clusteread.parse_fasta(ref_fasta)
    with open(os.path.join(path, "refs", ref + ".ref.fasta"), "w") as f:
        f.write(">" + id + "\n" + seq + "\n")

def read_hmmsearch_score(path):
    df = pd.read_csv(path,   delimiter=r"\s+", comment="#", usecols=[0,2,5], header=None)
    df.columns = ["sequence", "reference", "score"]
    df.sort_values(by=["score"], inplace=True, ascending=False)
    df["sequence"] = clean_column_ids(df, "sequence")
    df.set_index("sequence", inplace=True)
    return df

def find_highest_taxid(ncbi, fasta):
    id2taxid = taxtree.read_cluster_taxids(fasta)
    tax_tree = ncbi.get_topology(id2taxid.values(), intermediate_nodes=False)
    highest_tax = taxtree.get_highest_taxonomic_id(id2taxid.values(), tax_tree)
    for id in id2taxid.keys():
        if id2taxid[id] == highest_tax:
            return id
    
    assert False, "Highest taxid not found..."

def read_pairwise_identities(path):
    lines = []
    ids = {}

    with open(path, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            if line[0] == "#":
                continue
            seq1 = line.split()[0].split("|")[1]
            seq2 = line.split()[1].split("|")[1]
            id_percentage = float(line.split()[2])
            if not seq1 in ids.keys():
                ids[seq1] = []
            if not seq2 in ids.keys():
                ids[seq2] = []
            ids[seq1].append(id_percentage)
            ids[seq2].append(id_percentage)
    return ids

def find_highlow(ids):
    highlow = 0.0
    highlow_seqid = ""
    for seqid in ids.keys():
        m = min(ids[seqid])
        if m > highlow:
            highlow = m
            highlow_seqid = seqid
        if m == highlow:
            if sum(ids[seqid]) > sum(ids[highlow_seqid]):
                highlow_seqid = seqid

    return highlow_seqid

def find_max(ids):
    m = 0.0
    m_seqid = ""
    for seqid in ids.keys():
        s = sum(ids[seqid])
        if s > m:
            m = s
            m_seqid = seqid
    return m_seqid

def clustering(db, path):
    rmf([path + "/phmmer", path + "/refs"])
    clusters_path = glob.glob(os.path.join(path, "*.clusters"))[0]
    clusters, key_map = clusteread.read_clusters(db, clusters_path, 0, 100000, use_taxids=False)
    assert len(clusters) > 0, "No clusters found"
    fasta_files = sorted(glob.glob(os.path.join(path, "fasta", "*")))
    clean_files = sorted(glob.glob(os.path.join(path, "clean", "*")))
    files = zip(fasta_files, clean_files)
    for file_zip in files:
        with open(file_zip[0], "r") as f:
            oldref = f.readline().split(" ")[0][1:]
        reference = key_map[oldref].split("|")[1]
        set_refs(path, reference, file_zip)

# reqs: id/
def highlow(path):
    rmf([path + "/phmmer", path + "/refs"])
    id_files = sorted(glob.glob(os.path.join(path, "id", "*")))
    fasta_files = sorted(glob.glob(os.path.join(path, "fasta", "*")))
    clean_files = sorted(glob.glob(os.path.join(path, "clean", "*")))
    files = zip(fasta_files, clean_files, id_files)
    for file_zip in files:
        ids = read_pairwise_identities(file_zip[2])
        reference = find_highlow(ids)
        set_refs(path, reference, file_zip)

# reqs: id/
def identity(path):
    rmf([path + "/phmmer", path + "/refs"])
    id_files = sorted(glob.glob(os.path.join(path, "id", "*")))
    fasta_files = sorted(glob.glob(os.path.join(path, "fasta", "*")))
    clean_files = sorted(glob.glob(os.path.join(path, "clean", "*")))
    files = zip(fasta_files, clean_files, id_files)
    for file_zip in files:
        ids = read_pairwise_identities(file_zip[2])
        reference = find_max(ids)
        set_refs(path, reference, file_zip)

# reqs: hmmsearch/
def similarity(path):
    rmf([path + "/phmmer", path + "/refs"])
    hmms_files = sorted(glob.glob(os.path.join(path, "hmmsearch", "*")))
    fasta_files = sorted(glob.glob(os.path.join(path, "fasta", "*")))
    clean_files = sorted(glob.glob(os.path.join(path, "clean", "*")))
    files = zip(fasta_files, clean_files, hmms_files)
    for file_zip in files:
        hmms = read_hmmsearch_score(file_zip[2])
        reference = hmms.iloc[0]["reference"]
        set_refs(path, reference, file_zip)


def taxonomy(path):
    rmf([path + "/phmmer", path + "/refs"])
    fasta_files = sorted(glob.glob(path + "fasta", "*"))
    clean_files = sorted(glob.glob(path + "clean", "*"))
    files = zip(fasta_files, clean_files)
    ncbi = NCBITaxa()
    for file_zip in files:
        reference = find_highest_taxid(ncbi, file_zip[0])
        set_refs(path, reference.split("|")[1], file_zip)



if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("db", type=str, help="Database in fasta-format")
    parser.add_argument("path", type=str, help="Path to the cluster group directory")
    parser.add_argument("--clustering", action="store_true", help="First node in multi-step/mcl")
    parser.add_argument("--identity", action="store_true", help="Highest mean pair-wise identity score")
    parser.add_argument("--highlow", action="store_true", help="Highest lowest pair-wise identity score")
    parser.add_argument("--similarity", action="store_true", help="Highest hmmsearch score (i.e. most similar sequence to the MSA of the cluster)")
    parser.add_argument("--taxonomy", action="store_true", help="Highest node in taxonomic tree")
    args = parser.parse_args()
    if args.clustering:
        clustering(args.db, args.path)
    elif args.identity:
        identity(args.path)
    elif args.similarity:
        similarity(args.path)
    elif args.taxonomy:
        taxonomy(args.path)
    elif args.highlow:
        highlow(args.path)
    else:
        print("Reference criterion not specified. Use '-h' to list available criteria.")
