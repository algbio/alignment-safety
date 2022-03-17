def parse_fasta(protein_fasta):
    id = protein_fasta.split()[0]
    sequence = ''.join(protein_fasta.split("\n")[1:])
    return (id, sequence)

ids = []
with open("./../data/pdbs.txt", "r") as f:
    ids = f.read().splitlines()

with open("./../data/uniprot_sprot.fasta", "r") as f:
    db_fasta = ("\n" + f.read()).split("\n>")[1:]
    for protein in db_fasta:
        id, seq = parse_fasta(protein)
        accession = id.split("|")[1]
        if accession in ids:
            print(f">{protein}")

print(ids)