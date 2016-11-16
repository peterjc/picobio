#!/usr/bin/env python
"""Python script to insert gap features into EMBL files for ENA submission.

Given a FASTA assembly with runs of Ns in the sequence fed through Prokka
for annotation and then gff3_to_embl to make an EMBL file, it will fail
ENA validation:

$ java -jar embl-api-validator-1.1.146.jar BXF1.embl
...
ERROR: Sequence contains a stretch of 'n' characters between base {0} and {1} that is not represented with a "gap" feature (stretches of n greater than {2} gives a warning, greater than {3} gives an error). (26 occurrences) (SequenceToGapFeatureBasesCheck-1) 
...

References:

Prokka: https://github.com/tseemann/prokka

GFF3 to ENA ready EMBL: https://github.com/sanger-pathogens/gff3toembl

ENA Validator: https://www.ebi.ac.uk/ena/software/flat-file-validator
and https://github.com/enasequence/sequencetools

This script: https://github.com/peterjc/picobio/tree/master/seq_manipulation
"""

import sys

if len(sys.argv) != 3:
   sys.stop("Expects two arguments: Input EMLB filename, output EMBL filename\n")
input_embl = sys.argv[1]
output_embl = sys.argv[2]

try:
    from Bio import SeqIO
except ImportEror:
   sys.stop("This script requires Biopython 1.69 or later")


def insert_gaps(record):
   seq = str(record.seq).upper()
   sys.stderr.write("Record %s (length %i bp) has %i N characters\n" %
                    (record.id, len(seq), seq.count("N")))
   gap = "N" * 100
   try:
      i = seq.find(gap)
   except InderError:
      sys.stderr.write("No long gaps in record %s (length %i bp)\n" %
                       (record.id, len(seq)))
      return record

   count = 1
   sys.stderr.write("Added %i gap features to record %s (length %i bp)\n" %
                    (count, record.id, len(seq)))
   return record


print("Adding gap features to %s making %s" % (input_embl, output_embl))

# This is a memory efficient generator expression, one record at a time:
fixed_records = (insert_gaps(r) for r in SeqIO.parse(input_embl, "embl"))
count = SeqIO.write(fixed_records, output_embl, "embl")

print("Done, %i records" % count)
