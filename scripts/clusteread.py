import sys
import argparse

clusters = {}


def print_usage():
    print("Print cluster:   clusterread.py p <clustered db> <cluster file> <cluster id>")
    print("Write cluster:   clusterread.py w <clustered db> <cluster file> <cluster id>")
    print("Print info:      clusterread.py i <clustered db> <cluster file>")

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

            clusters[names[1]].append(names[0])
        print("Cluster count: %d" % cluster_count)
    finally:
        if mcl:
            for key in list(clusters.keys()):
                new_key = clusters[key][0]
                clusters[new_key] = clusters.pop(key)
        f.close()

def print_info():
    s_len, s_name = 1000000,""
    b_len, b_name = 0, ""
    for key in clusters.keys():
        if s_len > len(clusters[key]):
            s_len = len(clusters[key])
            s_name = key
        if b_len < len(clusters[key]):
            b_len = len(clusters[key])
            b_name = key

    print(f"Total number of clusters:       {len(clusters.keys())}")
    print(f"Size of the smallest cluster:   {s_len}, {s_name}")
    print(f"Size of the biggest cluster:    {b_len},  {b_name}")

def separate_clusters(db_filename, cluster_key):
    f = open(db_filename, "r")
    out = open("./out/" + cluster_key + ".fasta", "w")
    out_clean = open("./out/" + cluster_key + ".clean.fasta", "w")
    try:
        proteins = f.read().split("\n>")
        # new line at the beginning of the file so file can be split with "\n>"
        out.write("\n")
        for protein in proteins:
            protein_data = protein.split("\n")
            protein_id = protein_data[0].split(" ")[0]
            seq = ''.join(protein_data[1:])
            if protein_id in clusters[cluster_key]:
                out_clean.write(">" + protein_id + "\n" + seq + "\n")
                out.write(">" + protein + "\n")
    finally:
        f.close()

    out.close()
    out_clean.close()

        
def main():
    print(sys.argv)
    if len(sys.argv) < 4:
        print_usage()
        return
    read_clusters(sys.argv[3])
    if len(clusters) < 1:
        print("No clusters read...")
        return

    if sys.argv[1] == "i" and len(sys.argv) == 4:
        print_info()
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
            separate_clusters(sys.argv[2], sys.argv[4])
            print("Fasta file written succesfully...")
        else:
            print("No cluster found with cluster id: " + sys.argv[4])
    else:
        print_usage()


if __name__ == '__main__':
    main()
