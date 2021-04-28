import sys

from ete3 import Tree

cluster_name = ""

def print_stats(filename):
    t = Tree(filename)
    degree_of_branching = 0
    print(f"Stats of cluster: {cluster_name}")
    print(f"Number of sequences: {len(t.get_leaves())}")
    for node in t.traverse("postorder"):
        degree_of_branching += len(node.children) > 1
    print(f"Degree of branching: {degree_of_branching}")



if __name__ == '__main__':
    if len(sys.argv) == 2:
        path = sys.argv[1]
        # remove ./ from path
        if(path[:2] == "./"):
            path = path[2:]
        cluster_name = path.split("/")[1].split(".")[0]
        print_stats(path)