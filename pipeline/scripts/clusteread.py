import sys
import os
import argparse
# from ete3 import NCBITaxa
import random
from taxtree import get_highest_taxonomic_id, read_cluster_taxids


def read_clusters(db_file, filename, min_size, max_size):
    mcl = ".mcl" in filename
    mmseqs = ".mmseqs" in filename
    key_index = 0 if mmseqs else 1
    val_index = 1 if mmseqs else 0
    clusters = {}
    key_map = {}
    with open(filename, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break

            names = line.split()
            if names[key_index] == "-1":
                continue
            # if cluster does not exist yet
            if not names[key_index] in clusters:
                clusters[names[key_index]] = []
            
            clusters[names[key_index]].append(names[val_index])
    print(f"Total number of clusters in DB: {len(clusters.keys())}")

    # Delete clusters that dont fit min-max criteria
    for key in list(clusters.keys()):
        if min_size > len(clusters[key]) or len(clusters[key]) > max_size:
            del clusters[key]
    if mcl:
        for key in list(clusters.keys()):
                new_key = clusters[key][0]
                clusters[new_key] = clusters.pop(key)

    for key in clusters.keys():
        for prot in clusters[key]:
            key_map[prot] = key
    
    return (clusters, key_map)


def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

# Separates all clusters
def separate_clusters(clusters, key_map, db_filename, clustering_path, min_size, max_size, n=-1):
    mkdir(os.path.join(clustering_path, "fasta"))
    mkdir(os.path.join(clustering_path, "clean"))
    mkdir(os.path.join(clustering_path, "refs"))

    included = clusters.keys()
    if n > 0:
        included = random.sample(clusters.keys(), min(len(clusters.keys()), n))
            
    with open(os.path.join(clustering_path, "info.txt"), "w") as f:
        f.write(f"Database: {db_filename}\n")
        f.write(f"Clustering parameters: {clustering_path}\n")
        f.write(get_info(clusters))
        f.write(f"Cluster size range treshold: {min_size}-{max_size}\n")

    c = 0
    db_fasta = ""
    with open(db_filename, "r") as f:
        db_fasta = ("\n" + f.read()).split("\n>")[1:]
    
    cluster_num = {}
    ii = 1
    for key in clusters.keys():
        if key in included:
            cluster_num[key.split("|")[1]] = str(ii)
            ii += 1

    print("Writing reference sequences to fasta-files...")
    agh = []
    for protein_fasta in db_fasta:
        id, sequence = parse_fasta(protein_fasta)
        if id in included and id == key_map[id]:
            cleaned = id.split("|")[1]
            c += 1
            agh.append(id)
            # sys.stdout.write("\r%d%%" % int(c * 100.0 / len(included)))
            with open(os.path.join(clustering_path, "fasta", cluster_num[cleaned] + f"_{cleaned}" + ".fasta"), "w") as out:
                out.write(">" + protein_fasta + "\n")
            with open(os.path.join(clustering_path, "clean", cluster_num[cleaned] + f"_{cleaned}" + ".clean.fasta"), "w") as out:
                out.write(">" + id + "\n" + sequence + "\n")
            with open(os.path.join(clustering_path, "refs", cluster_num[cleaned] + f"_{cleaned}" + ".ref.fasta"), "w") as out:
                out.write(">" + id + "\n" + sequence + "\n")

    print("\nSeparating clusters to fasta-files...")
    for protein_fasta in db_fasta:
        if not protein_fasta:
            continue
        id, sequence = parse_fasta(protein_fasta)
        if id in key_map.keys() and key_map[id] in included and id != key_map[id]:
            # if(sequence.count("X") > 0):
            #     continue
            cluster_id = key_map[id]
            cleaned = cluster_id.split("|")[1]
            c += 1
            if not key_map[id] in agh:
                print(f"cluster: {id} not found")
            # sys.stdout.write("\r%d%%" % int(c * 100.0 / len(key_map.keys())))
            with open(os.path.join(clustering_path, "fasta", cluster_num[cleaned] + f"_{cleaned}" + ".fasta"), "a") as out:
                out.write(">" + protein_fasta + "\n")
            with open(os.path.join(clustering_path, "clean", cluster_num[cleaned] + f"_{cleaned}" + ".clean.fasta"), "a") as out:
                out.write(">" + id + "\n" + sequence + "\n")
    print("")
    with open(os.path.join(clustering_path, "info.txt"), "a") as f:
        f.write(f"Total number of sequences: {c}\n")
        f.write(f"Total number of clusters: {len(included)}\n")

            
def parse_fasta(protein_fasta):
    id = protein_fasta.split()[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)

# Separates one cluster
def separate_cluster(clusters, db_filename, cluster_key):
    db_fasta = ""
    with open(db_filename, "r") as f:
        db_fasta = ("\n" + f.read()).split("\n>")[1:]

    with open(os.path.join("out", cluster_key + ".fasta"), "w") as out:
        for protein in db_fasta:
            protein_id, sequence = parse_fasta(protein)
            if protein_id in clusters[cluster_key]:
                out.write(">" + protein + "\n")
    with open(os.path.join("out", cluster_key + ".clean.fasta"), "w") as out:
        for protein in db_fasta:
            protein_id, sequence = parse_fasta(protein)
            if protein_id in clusters[cluster_key]:
                out.write(">" + protein + "\n")

def get_info(clusters):
    info = ""
    s_len, s_name = 1000000,""
    b_len, b_name = 0, ""
    for key in clusters.keys():
        if s_len > len(clusters[key]):
            s_len = len(clusters[key])
            s_name = key
        if b_len < len(clusters[key]):
            b_len = len(clusters[key])
            b_name = key

    # bar diagram
    # 1-10, 11-100, 101-300, 301-500, 501-1000, 1000+
    ranges = [10, 100, 300, 500, 1000]
    bins = [0,0,0,0,0,0]
    for key in clusters.keys():
        l = len(clusters[key])
        if 0 < l <= ranges[0]:
            bins[0] += 1
        elif ranges[0] < l <= ranges[1]:
            bins[1] += 1
        elif ranges[1] < l <= ranges[2]:
            bins[2] += 1
        elif ranges[2] < l <= ranges[3]:
            bins[3] += 1
        elif ranges[3] < l <= ranges[4]:
            bins[4] += 1
        elif ranges[4] < l:
            bins[5] += 1

    info += f"Total number of clusters:       {len(clusters.keys())}\n"
    info += f"Size of the smallest cluster:   {s_len}, {s_name}\n"
    info += f"Size of the biggest cluster:    {b_len},  {b_name}\n"
    info += "Distribution of cluster sizes:\n"
    info += "Size range  : Cluster count\n"
    for i in range(6):
        if i == 0:
            info += f"{0:>4} - {ranges[i]:>4} : {bins[i]:6}\n"
        elif i == 5:
            info += f"{ranges[i-1]+1:>4} +      : {bins[i]:6}\n"
        else:
            info += f"{ranges[i-1]+1:>4} - {ranges[i]:>4} : {bins[i]:6}\n"
    return info

def main(args):
    clusters, key_map = read_clusters(args.db, args.clusters, args.min, args.max)
    if len(clusters) < 1:
        print("No clusters read...")
        return

    if args.action == "a":
        separate_clusters(clusters, key_map, args.db, "/".join(args.clusters.split("/")[:-1]), args.min, args.max, args.n)

    elif args.action == "i":
        print(get_info(clusters))

    elif args.action == "p":
        if not args.id in clusters:
            print("Cluster not found...")
            return
        cluster_size = 0
        for protein in clusters[args.id]:
            print(protein)
            cluster_size += 1
        print("Cluster size: %d" % cluster_size)

    elif args.action == "w":
        print("Writing to file...")
        if args.id in clusters:
            separate_cluster(clusters, args.db, args.id)
            print("Fasta file written succesfully...")
        else:
            print("No cluster found with cluster id: " + args.id)


def check_args(parser, args):
    if args.action == "a":
        return True
    elif args.action == "i":
        return True
    elif args.action == "p":
        if not "id" in vars(args):
            parser.error("Action \"p\" requires --id")
        return "id" in vars(args)
    elif args.action == "w":
        if not "id" in vars(args):
            parser.error("Action \"w\" requires --id")
        return "id" in vars(args)
    else:
        return False

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("action", type=str, help="a/i/p/w : separate all/print info/print cluster/write cluster")
    parser.add_argument("db", type=str, help="database in fasta-format")
    parser.add_argument("clusters", type=str, help="cluster file")
    parser.add_argument("--min", type=int, default=10, help="minimum size of the cluster {10}")
    parser.add_argument("--max", type=int, default=10000, help="maximum size of the cluster {1000}")
    parser.add_argument("--n", type=int, default=-1, help="If specified, n-number of random clusters will be selected")
    parser.add_argument("--id", type=str, help="cluster id to separate")
    args = parser.parse_args()
    if check_args(parser, args):
        main(args)
