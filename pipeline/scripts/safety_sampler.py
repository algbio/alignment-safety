import argparse
import random

def read_safety_windows(file):
    with open(file, "r") as f:
        data = {}
        seqID = {} # Sequence of IDs per cluster
        seq = {} # Sequence per cluster
        references = {} # Reference sequence for each cluster
        referencesID = {} # Reference ID sequence for each cluster
        nS = {} # Number of sequence in each cluster
        intervals = {} # Interval for each safety window in each sequence/reference in each cluster
        intervalsAligned = {} # Intervals for each safety in each sequence in each cluster
        windowC = {} # Number of windows per sequences per cluster
        
        fn = open(file,'r') # read the file
            
        line = fn.readline().replace('"', '').split() # reference line 
        referencesID = line[0] # reference ID sequence in that cluster
        references = line[1] # reference sequence in that cluster
        nS = int(fn.readline().splitlines()[0]) # number of other sequence in that cluster
        for s in range(1,nS):
            line = fn.readline().replace('"', '').split()
            seqID[s] = str(line[0])
            seq[s] = str(line[1])
            intervals[s] = list()
            intervalsAligned[s] = list()
            for k in range(0,int(line[2])):
                line = fn.readline().split()
                intervals[s].append((int(line[0]),int(line[1])))
                intervalsAligned[s].append((int(line[2]),int(line[3])))
        
        data['reference'] = references
        data['sequences per cluster'] = nS
        data['safety windows intervals'] = intervals
        data['sequences windows intervals'] = intervalsAligned
        data['sequences'] = seq
        
        return data
    assert False, "Reading file not succesful"

def generate_random_sequences(data, output_file):
    with open(output_file, "w") as out:
        for seq_i in data["sequences"].keys():
            out.write(f">rs{seq_i}\n")
            random_seq = sample_from(data["sequences"][seq_i], len(data["sequences"][seq_i]))
            for window in data["sequences windows intervals"][seq_i]:
                start = window[0]
                end = window[1]
                random_seq = random_seq[:start] + data["sequences"][seq_i][start:end+1] + random_seq[end+1:]
            out.write(random_seq + "\n")

def sample_from(seq, l):
    seqlen = len(seq)-1
    out = ""
    for i in range(l):
        r = random.randint(0, seqlen)
        out += seq[r]
    return out

def main(args):
    data = read_safety_windows(args.safety_file)
    print(data)
    generate_random_sequences(data, args.output_file)

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("safety_file", type=str, help="Emerald output file")
    parser.add_argument("output_file", type=str, help="Where to write output (in fasta-format)")
    parser.add_argument("--inverse", action="store_true", help="To only randomly sample inside of safety intervals")
    args = parser.parse_args()
    main(args)