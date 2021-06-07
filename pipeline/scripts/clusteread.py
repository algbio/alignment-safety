import os
import sys
import argparse

clusters = {}
key_map = {}


def belongs_to_cluster(protein_id):
    return key_map[protein_id]

def print_usage():
    print("Print cluster:   clusterread.py p <original db> <cluster file> <cluster id>")
    print("Print cluster:   clusterread.py a <original db> <cluster file> <min size> <max size>")
    print("Write cluster:   clusterread.py w <original db> <cluster file> <cluster id>")
    print("Print info:      clusterread.py i <original db> <cluster file>")

def read_clusters(filename):
    # check which file format clusters are
    mcl = filename.split(".")[-1] == "mcl"
    f = open(filename, "r")
    try:
        cluster_count = 0
        lines = f.readlines()
        for line in lines:
            names = line.split()
            # if cluster does not exist yet
            if not names[1] in clusters:
                clusters[names[1]] = []
                cluster_count += 1

            key_map[names[0]] = names[1]
            clusters[names[1]].append(names[0])
        print("Cluster count: %d" % cluster_count)
    finally:
        if mcl:
            for key in list(clusters.keys()):
                new_key = clusters[key][0]
                clusters[new_key] = clusters.pop(key)
        f.close()

def get_info():
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

def separate_clusters(db_filename, clustering_path, min_size, max_size):
    print(clustering_path)
    if not os.path.exists(clustering_path + "/fasta"):
        os.makedirs(clustering_path + "/fasta")
    if not os.path.exists(clustering_path + "/clean"):
        os.makedirs(clustering_path + "/clean")
    with open(clustering_path + "/info.txt", "w") as f:
        f.write(f"Database: {db_filename}\n")
        f.write(f"Clustering parameters: {clustering_path}\n")
        f.write(get_info())
        f.write(f"Cluster size range treshold: {min_size}-{max_size}")
    with open(db_filename, "r") as f:
        db_fasta = ("\n" + f.read()).split("\n>")
        i = 0
        for protein_fasta in db_fasta:
            if i % 2000 == 0:
                print(f"progress: {i * 100.0 / len(db_fasta):.1f}%")
            if not protein_fasta:
                continue
            id, sequence = parse_fasta(protein_fasta)
            cluster_id = belongs_to_cluster(id)
            if(min_size <= len(clusters[cluster_id]) <= max_size):
                with open(clustering_path + "/fasta/" + cluster_id + ".fasta", "a") as out:
                    out.write(">" + protein_fasta + "\n")
                with open(clustering_path + "/clean/" + cluster_id + ".clean.fasta", "a") as out:
                    out.write(">" + id + "\n" + sequence + "\n")
            i += 1

            
def parse_fasta(protein_fasta):
    id = protein_fasta.split(" ")[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)

def separate_cluster(db_filename, cluster_key):
    f = open(db_filename, "r")
    out = open("./out/" + cluster_key + ".fasta", "w")
    out_clean = open("./out/" + cluster_key + ".clean.fasta", "w")
    try:
        proteins = f.read().split("\n>")
        # new line at the beginning of the file so file can be split with "\n>"
        out.write("\n")
        for protein in proteins:
            protein_id, sequence = parse_fasta(protein)
            if protein_id in clusters[cluster_key]:
                out_clean.write(">" + protein_id + "\n" + sequence + "\n")
                out.write(">" + protein + "\n")
    finally:
        f.close()

    out.close()
    out_clean.close()

        
def main():
    if len(sys.argv) < 4:
        print_usage()
        return

    read_clusters(sys.argv[3])

    if len(clusters) < 1:
        print("No clusters read...")
        return

    if sys.argv[1] == "a" and len(sys.argv) == 6:
        separate_clusters(sys.argv[2], sys.argv[3] + ".clusters", int(sys.argv[4]), int(sys.argv[5]))
        return

    elif sys.argv[1] == "i" and len(sys.argv) == 4:
        print(get_info())
        return

    elif sys.argv[1] == "p" and len(sys.argv) == 5:
        if not sys.argv[3] in clusters:
            print("Cluster not found...")
            return
        
        cluster_size = 0
        for protein in clusters[sys.argv[4]]:
            print(protein)
            cluster_size += 1
        print("Cluster size: %d" % cluster_size)

    elif sys.argv[1] == "w" and len(sys.argv) == 5:
        print("Writing to file...")
        if sys.argv[4] in clusters:
            separate_cluster(sys.argv[2], sys.argv[4])
            print("Fasta file written succesfully...")
        else:
            print("No cluster found with cluster id: " + sys.argv[4])
    else:
        print_usage()


if __name__ == '__main__':
    main()
