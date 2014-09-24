#!/usr/bin/env python
"""Python script for building Biopython SeqIO SQLite index files.

Intended for use as part of a larger pipeline, e.g. with make to
ensure every input FASTQ file has been indexed.
"""

import os
import sys
from optparse import OptionParser

VERSION = "0.0.1"

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

try:
    from Bio.SeqIO import index_db
except ImportError:
    sys_exit("Missing Biopython (or too old to provide Bio.SeqIO.index_db(...)")


def main():
    usage = """usage: %prog [options] [sequence filenames]

    By default (if not using the -i or --index argument), one index
    will be created for each sequence file with ".idx" appended.
    """
    parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog "+VERSION)
    parser.add_option("-f", "--format", dest="format",
                      type="string", metavar="FORMAT",
                      help="""Sequence format supported by Biopython's SeqIO.

                      Common examples would be 'fasta' or 'fastq' (required).
                      """)
    parser.add_option("-i", "--index", dest="index_filename",
                      type="string", metavar="FILE",
                      help="Created one combined index using this filename.")
    (options, args) = parser.parse_args()

    if not args:
        sys_exit("No sequence filenames provided")
    if not options.format:
        sys_exit("No sequence format specified")
    format = options.format.lower()

    for filename in args:
        if not os.path.isfile(filename):
            sys_exit("Missing %s" % filename)
    
    if options.index_filename:
        # One shared index for all sequence files
        idx_filename = options.index_filename
        if not os.path.isfile(idx_filename):
            print("%s - indexing %i files..." % len(filenames))
        d = index_db(idx_filename, filenames, format)
        print("%s - OK, %i records in %i files" % (idx_filename, len(d), len(filenames)))
    else:
        # One index per sequence file
        for filename in args:
            idx_filename = filename + ".idx"
            if not os.path.isfile(idx_filename):
                print("%s - indexing..." % filename)
            d = index_db(idx_filename, filename, format)
            print("%s - OK, %i records" % (idx_filename, len(d)))


if __name__ == "__main__":
    main()
