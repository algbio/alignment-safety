import sys
from os.path import expanduser
import subprocess as sp
from ete3 import Tree, NCBITaxa
ncbi = NCBITaxa()
from temp import get_time_ms

# Downloads/Updates taxonomy database
# force_update = True
# if(force_update):
#     ncbi.update_taxonomy_database("/wrk/users/gyntartu/pipeline/data/ncbi")
# else:
#     try:
#         file = str(expanduser("~")) + "/.etetoolkit/taxa.sqlite"
#         f = open(file)
#         f.close()
#     except IOError:
#         ncbi.update_taxonomy_database()
    

ids_to_taxid = {}
cluster_name = ""

def get_highest_taxonomic_id(taxids, tax_tree):
    for node in tax_tree.traverse("levelorder"):
        if node.taxid in taxids:
            return node.taxid
    return taxids[0]

def subtree():
    ids = list(ids_to_taxid.values())
    tax_tree = ncbi.get_topology(ids, intermediate_nodes=True)
    for protein_id in ids_to_taxid.keys():
        tax_id = ids_to_taxid[protein_id]
        parent = tax_tree.search_nodes(taxid=tax_id)[0]
        parent.add_child(name=protein_id)

    print_stats(tax_tree)
    print(tax_tree.get_ascii(attributes=["taxid", "name"]))
    print(get_highest_taxonomic_id(ids, tax_tree))
    # tax_tree.write(features=[], format=2, outfile="./out/" + cluster_name + ".tax.tree")

def print_stats(tree):
    rank_counts = {}
    degree_of_branching = 0
    print(f"Stats of cluster: {cluster_name}")
    print(f"Number of sequences: {len(ids_to_taxid.keys())}")
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
    

def read_cluster_taxids(filename):
    tax_dict = {}
    proteins = ""
    with open(filename, "r") as f:
        proteins = ("\n" + f.read()).split("\n>")[1:]

    for protein in proteins:
        protein_data = protein.split("\n")[0]
        protein_id = protein_data.split(" ")[0]
        slc = protein_data.find("OX=")
        tax_id = int(protein_data[slc+3:].split(" ")[0])
        try:
            ncbi.get_lineage(tax_id)
            tax_dict[protein_id] = tax_id
        except:
            print(f"Unknown taxid: {tax_id}")

    return tax_dict
    

if __name__ == '__main__':
    if len(sys.argv) > 1:
        path = sys.argv[1]
        # remove ./ from path
        if(path[:2] == "./"):
            path = path[2:]
        cluster_name = path.split("/")[1].split(".")[0]
        ids_to_taxid = read_cluster_taxids(path)
        subtree()