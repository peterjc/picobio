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

blast_file = assembly_fasta + ".blast.tsv"

if not os.path.isfile(assembly_fasta):
    stop_err("Assemlby FASTA file not found: %r" % assembly_fasta)

if not os.path.isfile(reference_fasta):
    stop_err("Reference FASTA file not found: %r" % reference_fasta)

#TODO - Actual code...
def do_blast(query_fasta, db_fasta, blast_file):
    assert os.path.isfile(query_fasta)
    assert os.path.isfile(db_fasta)
    assert os.path.isfile(db_fasta + ".nhr")
    assert os.path.isfile(db_fasta + ".nin")
    assert os.path.isfile(db_fasta + ".nsq")
    cmd = NcbiblastnCommandline(query=query_fasta, db=db_fasta,
                                out=blast_file, outfmt=6,
                                evalue=1e-5)
    print cmd
    stdout, stderr = cmd()
    return

if not os.path.isfile(blast_file):
    do_blast(assembly_fasta, reference_fasta, blast_file)
