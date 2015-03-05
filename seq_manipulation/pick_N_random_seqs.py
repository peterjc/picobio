#!/usr/bin/env python
import random
import sys
from Bio import SeqIO
from Bio.SeqIO.FastaIO import SimpleFastaParser

"""Extract N randomly selected sequences from a FASTA file.

Takes three arguments: input FASTA filename, number of
sequences to pick out, and output FASTA filename. e.g.

$ python pick_N_random_seqs.py input.fasta 1000 output.fasta

If the input FASTA file has less than the requested count,
this will fail with an error.
"""

input_fasta, count, output_fasta = sys.argv[1:]
count = int(count)

with open(input_fasta) as handle:
    # Using as faster than SeqIO.parse(...)
    ids = [title.split(None, 1)[0] for title, seq in SimpleFastaParser(handle)]
print("Input FASTA file %s has %i sequences"
      % (input_fasta, len(ids)))
assert len(set(ids)) == len(ids), "You have duplicate identifiers"

#seqs = SeqIO.index(input_fasta, "fasta")
#print("Input FASTA file %s has %i sequences"
#      % (input_fasta, len(seqs)))
#assert count <= len(seqs)
#picked = set(random.sample(list(seqs), count))
#assert len(picked) == count
#del seqs

picked = set(random.sample(ids, count))

# This will preserve the input order, and do line wrapping
wanted = (r for r in SeqIO.parse(input_fasta, "fasta") if r.id in picked)
saved = SeqIO.write(wanted, output_fasta, "fasta")
assert saved == count

print("Saved %i randomly selected records from %s into %s"
      % (count, input_fasta, output_fasta))
