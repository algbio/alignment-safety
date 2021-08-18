#!/usr/bin/env python3

from Bio import SeqIO
        
contig_count = 0
with open(snakemake.input.assembly, 'r') as infile:
    for record in SeqIO.parse(infile, 'fasta'):
        contig_count += 1

with open(snakemake.output.report, 'w') as outfile:
    outfile.write('contig_count: {}\n'.format(contig_count))