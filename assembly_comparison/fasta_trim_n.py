#!/usr/bin/env python
"""Python script for trimming N bases from ends of sequences.
"""

import os
import sys
from optparse import OptionParser

usage = """Basic usage: ./fasta_trim_n.py < input.fasta > output.fasta

For more details, run with -h for the help.
"""

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

try:
    from Bio import SeqIO
except ImportError:
    stop_err("This script requires Biopython")

parser = OptionParser(usage=usage)
parser.add_option("-i", "--input", dest="input_filename",
                  help="Input sequence file (default is stdin)",
                  default=None,
                  metavar="FILE")
parser.add_option("-o", "--output", dest="output_filename",
                  help="Output sequence file (fefault is stdout)",
                  default=None,
                  metavar="FILE")
parser.add_option("-f", "--format", dest="sequence_format",
                  help='Sequence format (as named in Biopython SeqIO, default "fasta")',
                  default="fasta")
parser.add_option("-c", "--chars", dest="characters",
                  help='Characters to trim (default "Nn" covering upper and lower case)',
                  default="Nn",
                  metavar="FILE")
(options, args) = parser.parse_args()

chars = options.characters
format = options.sequence_format.lower()

sys.stderr.write("Removing %s characters from start/end of %s format file...\n" % (chars, format))

if options.input_filename:
    input_handle = open(options.input_filename)
else:
    input_handle = sys.stdin

if options.output_filename:
    output_handle = open(options.output_filename, "w")
else:
    output_handle = sys.stdout

chars = options.characters
format = options.sequence_format.lower()

def strip_seq(records):
    for record in records:
        #FASTQ etc will be a problem, must trim quality too!
        #old_len = len(record.seq)
        record.seq = record.seq.strip(chars)
        #TODO Minium length!
        #new_len = len(record.seq)
        #if new_len < old_len:
        #    sys.stderr.write("Trimmed %s from %i to %i\n" % (record.id, old_len, new_len))
        yield record

#Do the work,
count = SeqIO.write(strip_seq(SeqIO.parse(input_handle,format)), output_handle, format)

if options.input_filename:
    input_handle.close()
if options.output_filename:
    output_handle.close()

sys.stderr.write("Saved %i records\n" %count)
