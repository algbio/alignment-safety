import sys
from os.path import expanduser
import subprocess as sp
from ete3 import Tree, NCBITaxa

ncbi = NCBITaxa()

# Downloads/Updates taxonomy database
# Uncomment this line for first run only

force_update = False
if(force_update):
    ncbi.update_taxonomy_database()
else:
    try:
        file = str(expanduser("~")) + "/.etetoolkit/taxa.sqlite"
        f = open(file)
        f.close()
    except IOError:
        ncbi.update_taxonomy_database()
    

tax_dict = {}
cluster_name = ""

def subtree():
    ids = list(tax_dict.values())
    tax_tree = ncbi.get_topology(ids)
    for protein_id in tax_dict.keys():
        tax_id = tax_dict[protein_id]
        parent = tax_tree.search_nodes(taxid=tax_id)[0]
        parent.add_child(name=protein_id)

    print_stats(tax_tree)
    # print(tax_tree.get_ascii(attributes=["rank", "name"]))
    tax_tree.write(features=[], format=2, outfile="./out/" + cluster_name + ".tax.tree")

def print_stats(tree):
    rank_counts = {}
    degree_of_branching = 0
    print(f"Stats of cluster: {cluster_name}")
    print(f"Number of sequences: {len(tax_dict.keys())}")
    for node in tree.traverse("postorder"):
        if hasattr(node, "rank"):
            if not node.rank in rank_counts.keys():
                rank_counts[node.rank] = 0
            rank_counts[node.rank] += 1
        degree_of_branching += len(node.children) > 1
    print(f"Degree of branching: {degree_of_branching}")
    print("Number of distinct ranks:")
    for key in rank_counts.keys():
        print(f"Rank: {key} : {rank_counts[key]}")
    


def read_cluster_ids(filename):
    f = open(filename, "r")
    try:
        proteins = f.read().split("\n>")[1:]
        for protein in proteins:
            protein_data = protein.split("\n")[0]
            protein_id = protein_data.split(" ")[0]
            slc = protein_data.find("OX=")
            tax_id = int(protein_data[slc+3:].split(" ")[0])
            tax_dict[protein_id] = tax_id

    finally:
        f.close()
    

if __name__ == '__main__':
    if len(sys.argv) > 1:
        path = sys.argv[1]
        # remove ./ from path
        if(path[:2] == "./"):
            path = path[2:]
        cluster_name = path.split("/")[1].split(".")[0]
        read_cluster_ids(path)
        subtree()