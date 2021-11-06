import argparse

def read_safety_windows(file):
    with open(file, "r") as f:
        while True:
            line = f.readline()
            if not line:
                break
            
# read input files
def readSafety(data,file):
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
            intervals[s].append((line[0],line[1]))
            intervalsAligned[s].append((line[2],line[3]))
    
    data['reference'] = references
    data['sequences per cluster'] = nS
    data['safety windows intervals'] = intervals
    data['sequences windows intervals'] = intervalsAligned
    data['sequences'] = seq
    
    return data;

def main(args):
    pass

if __name__ == '__main__':
    parser = argparse.ArgumentParser()
    parser.add_argument("safety_file", type=str, help="Emerald output file")
    parser.add_argument("output_file", type=str, help="Where to write output (in fasta-format)")
    parser.add_argument("--inverse", action="store_true", help="To only randomly sample inside of safety intervals")
    args = parser.parse_args()
    main(args)