import os
import argparse
from posixpath import split
from unittest.util import safe_repr
import pandas as pd
import numpy as np

from Safe import *

def read_pdb(pdb):
    header = ""
    footer = ""
    aas = []
    with open(pdb, "r") as f:
        line = f.readline()
        while line.split()[0] != "ATOM":
            header += line
            line = f.readline()
        while line.split()[0] == "ATOM":
            cur_aa = line.split()[5]
            aa = ""
            while line.split()[0] == "ATOM" and cur_aa == line.split()[5]:
                aa += line
                line = f.readline()
            aas.append(aa)
        while line:
            footer += line
            line = f.readline()
    return (header, aas, footer)

def mkdir(path):
    if not os.path.exists(path):
        os.makedirs(path)

def write_aas(aas, out):
    with open(out, "w") as f:
        assert len(aas) > 3, "Should be more than 3 bases"
        for aa in aas:
            f.write(aa)

def split_pdb(ss, min, out, ref_aas, seq_aas):
    unsafe_path = os.path.join(out, ss.accession, "unsafe")
    safe_path = os.path.join(out, ss.accession, "safe")
    mkdir(unsafe_path)
    mkdir(safe_path)

    windows = list(zip(ss.raw_ref_windows, ss.raw_seq_windows))
    c_safe = 0
    c_unsafe = 0
    for (i, ((rstart, rend), (sstart, send))) in enumerate(windows):
        if i > 0:
            ((_, prend), (_, psend)) = windows[i-1]
            if rstart - prend > min and  sstart - psend > min:
                write_aas(ref_aas[prend:rstart], os.path.join(unsafe_path, f"ref_{c_unsafe}.pdb"))
                write_aas(seq_aas[psend:sstart], os.path.join(unsafe_path, f"seq_{c_unsafe}.pdb"))
                c_unsafe += 1

        if rend - rstart > min and send - sstart > min:
            write_aas(ref_aas[rstart:rend+1], os.path.join(safe_path, f"ref_{c_safe}.pdb"))
            write_aas(seq_aas[sstart:send+1], os.path.join(safe_path, f"seq_{c_safe}.pdb"))
            c_safe += 1


def main(args):
    cluster = Cluster(args.safe)
    _, ref_aas, _ = read_pdb(os.path.join(args.pdb,f"{cluster.ref_accession}.pdb"))
    for ss in cluster.sequences:
        _, seq_aas, _ = read_pdb(os.path.join(args.pdb,f"{ss.accession}.pdb"))
        split_pdb(ss, args.min, args.out, ref_aas, seq_aas)
    


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("pdb", type=str, help="PDB directory")
    parser.add_argument("safe", type=str, help="Safety file")
    parser.add_argument("out", type=str, help="Output directory")
    # parser.add_argument("--th", type=float, default=0.5, help="Threshold of Safeness")
    parser.add_argument("--min", type=int, default=10, help="Minimum number of amino acids to cut.")
    args = parser.parse_args()
    main(args)