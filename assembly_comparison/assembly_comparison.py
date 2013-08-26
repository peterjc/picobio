#!/usr/bin/env python
"""Python script for assembly comparison.
"""

import os
import sys
from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

#TODO - Use argparse if API becomes non-trivial.
if len(sys.argv) != 3:
    stop_err("Usage: do_comparison.py assembly.fasta reference.fasta")
assembly_fasta, reference_fasta = sys.argv[1:]

if not os.path.isfile(assembly_fasta):
    stop_err("Assemlby FASTA file not found: %r" % assembly_fasta)

if not os.path.isfile(reference_fasta):
    stop_err("Reference FASTA file not found: %r" % reference_fasta)

#TODO - Actual code...
