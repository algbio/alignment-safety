import glob
import os
import argparse

import clusteread

def validate_fasta(path):
    data = {}
    with open(path, "r") as f:
        db_fasta = ("\n" + f.read()).split("\n>")[1:]
        for protein in db_fasta:
            protein_id, sequence = clusteread.parse_fasta(protein)
            clean_id = protein_id.split("|")[1]
            assert not clean_id in data.keys(), f"Duplicate protein: {path} {clean_id}"
            data[clean_id] = protein
    return len(data.keys())

def validate_group(path):
    fasta_files = sorted(glob.glob(os.path.join(path, "fasta", "*.fasta")))
    clean_files = sorted(glob.glob(os.path.join(path, "clean", "*.clean.fasta")))
    ref_files = sorted(glob.glob(os.path.join(path, "refs", "*.ref.fasta")))
    assert (len(fasta_files) == len(clean_files) and len(fasta_files) == len(ref_files)), "Wrong amount of files"
    files = zip(fasta_files, clean_files, ref_files)
    fasta_count = 0
    clean_count = 0
    for filezip in files:
        fasta_count += validate_fasta(filezip[0])
        clean_count += validate_fasta(filezip[1])
    assert fasta_count == clean_count, "Wrong amount of sequences"
    print(f"Validation OK! Clusters: {len(fasta_files)}, Seqs: {fasta_count}")


if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("path", type=str, help="Group path")
    args = parser.parse_args()
    validate_group(args.path)