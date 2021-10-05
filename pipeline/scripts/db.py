def parse_fasta(protein_fasta):
    id = protein_fasta.split()[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)

db_fasta = ""
with open("data/uniprot_sprot.fasta", "r") as f:
    db_fasta = ("\n" + f.read()).split("\n>")[1:]

counts = {}
lens = {}
for protein in db_fasta:
    protein_id, sequence = parse_fasta(protein)
    if not sequence[0] in counts:
        lens[sequence[0]] = 0
        counts[sequence[0]] = 0
    counts[sequence[0]] += 1
    lens[sequence[0]] += len(sequence)
    
print(f"Total sequences: {len(db_fasta)}")
l = len(db_fasta)
for aa in counts.keys():
    print(f"{aa} : {counts[aa]} ({counts[aa]*100.0/l:.1f}%)", end=" ")
    print(f"Avg. seq length: {lens[aa] / counts[aa]:.1f}")