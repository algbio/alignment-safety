import sys
from os.path import expanduser
import subprocess as sp
from ete3 import Tree, NCBITaxa
import time
ncbi = NCBITaxa()



def get_highest_taxonomic_sequence(ids):
    tax_tree = ncbi.get_topology(ids)
    for node in tax_tree.traverse("levelorder"):
        print(node)


def get_time_ms():
    return round(time.time() * 1000)