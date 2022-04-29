import argparse
import random

from Safe import Cluster, merge_windows, read_safe_sequences


def generate_random_sequences(cluster, output_file, fill_x):
    with open(output_file, "w") as out:
        out.write(f">{cluster.ref_accession}\n")
        out.write(cluster.ref_seq + "\n")
        for seq in cluster.sequences:
            out.write(f">{seq.accession}\n")

            if fill_x:
                random_seq = "X" * len(seq.sequence)
            else:
                random_seq = sample_from(seq.sequence, len(seq.sequence))
                
            for window in merge_windows(seq.raw_seq_windows):
                start = window[0]
                end = window[1]
                random_seq = random_seq[:start] + seq.sequence[start:end+1] + random_seq[end+1:]

            out.write(random_seq + "\n")

def sample_from(seq, l):
    seqlen = len(seq)-1
    out = ""
    for i in range(l):
        r = random.randint(0, seqlen)
        out += seq[r]
    return out

def main(args):
    cluster = Cluster(args.safety_file)
    generate_random_sequences(cluster, args.output_file, args.x)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("safety_file", type=str, help="Emerald output file")
    parser.add_argument("output_file", type=str, help="Where to write output (in fasta-format)")
    # parser.add_argument("--inverse", action="store_true", help="To only randomly sample inside of safety intervals")
    parser.add_argument("-x", action="store_true", help="Fill safety intervals with 'X'")
    args = parser.parse_args()
    main(args)