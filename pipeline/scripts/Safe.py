import numpy as np

def parse_fasta(protein_fasta):
    id = protein_fasta.split("\n")[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)

def merge_windows(wins):
    a = [wins[0]]
    for i in range(len(wins)-1):
        prev = a[len(a)-1]
        next = wins[i+1]
        if prev[1] + 1 >= next[0]:
            a[len(a)-1] = (a[len(a)-1][0], next[1])
        else:
            a.append(next)
    return a

class SafeSequence():
    def __init__(self, accession, sequence, raw_ref_windows, raw_seq_windows):
        self.accession = accession
        self.sequence = sequence
        self.raw_ref_windows = raw_ref_windows
        self.raw_seq_windows = raw_seq_windows
        self.bases = np.zeros(len(sequence), dtype=bool)
        self.init()
        self.safe_coverage = np.sum(self.bases) * 1.0 / len(sequence)

    def init(self):
        for window in merge_windows(self.raw_seq_windows):
            start = window[0]
            end = window[1]
            for i in range(start, end):
                self.bases[i] = True

    def is_safe(self, i):
        return self.bases[i]

def read_safe_sequences(safety_path):
    safe_sequences = []
    raw = ""
    with open(safety_path, "r") as f:
        raw = ("\n" + f.read()).split("\n>")[1:]
    ref_accession = raw[0].split("\n")[0].split(" ")[0].strip().split("|")[1]
    ref = raw[0].split("\n")[1].strip()
    
    for data in raw[1:]:
        splitted = data.split("\n")
        if ":" in splitted[0].split(" ")[0]:
            accession = splitted[0].split(" ")[0].split(":")[1]
        else:
            accession = splitted[0].split(" ")[0].split("|")[1]
        seq = splitted[1].strip()
        rwindows = []
        swindows = []
        for line in splitted[3:]:
            ends = line.split()
            if len(ends) != 4:
                continue
            rwindows.append((int(ends[0]), int(ends[1])))
            swindows.append((int(ends[2]), int(ends[3])))
        safe_sequences.append(SafeSequence(accession, seq, rwindows, swindows))
        
    return (ref, ref_accession, safe_sequences)

class Cluster():
    def __init__(self, safe_path):
        self.ref_seq, self.ref_accession, self.sequences = read_safe_sequences(safe_path)
        self.mapping = {}
        for seq in self.sequences:
            self.mapping[seq.accession] = seq

    def get(self, accession):
        return self.mapping[accession]

            