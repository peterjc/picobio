#!/usr/bin/env python
"""Python script for building Biopython SeqIO SQLite index files.

Intended for use as part of a larger pipeline, e.g. with make to
ensure every input FASTQ file has been indexed.

History:

* v0.0.1 - Original
* v0.0.2 - Option to force re-creation of indexes
         - Check indexes claiming to have zero records
* v0.0.3 - Actually print the usage text I wrote.
"""

import os
import sys
from optparse import OptionParser

VERSION = "0.0.3"


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

try:
    from Bio import SeqIO
except ImportError:
    sys_exit("Missing Biopython")

try:
    from Bio.SeqIO import index_db
except ImportError:
    sys_exit("Biopython too old to provide Bio.SeqIO.index_db(...)?")


def at_least_one_record(filenames, format):
    for f in filenames:
        for r in SeqIO.parse(f, format):
            return True
    # Really were no records!
    return False


def main():
    usage = """usage: %prog [options] [sequence filenames]

    By default (if not using the -i or --index argument), one index
    will be created for each sequence file with ".idx" appended.
    """
    parser = OptionParser(usage=usage,
                          version="%prog " + VERSION)
    parser.add_option("-f", "--format", dest="format",
                      type="string", metavar="FORMAT",
                      help="""Sequence format supported by Biopython's SeqIO.

                      Common examples would be 'fasta' or 'fastq' (required).
                      """)
    parser.add_option("-i", "--index", dest="index_filename",
                      type="string", metavar="FILE",
                      help="Created one combined index using this filename.")
    parser.add_option("-r", "--reindex", dest="reindex",
                      action="store_true",
                      help="Delete pre-existing indexes and rebuild them.")
    (options, filenames) = parser.parse_args()

    if not filenames:
        sys_exit("No sequence filenames provided")
    if not options.format:
        sys_exit("No sequence format specified")
    format = options.format.lower()

    for filename in filenames:
        if not os.path.isfile(filename):
            sys_exit("Missing %s" % filename)

    if options.index_filename:
        # One shared index for all sequence files
        idx_filename = options.index_filename
        if options.reindex and os.path.isfile(idx_filename):
            print("%s - re-indexing %i files.." % len(filenames))
            os.remove(idx_filename)
        elif not os.path.isfile(idx_filename):
            print("%s - indexing %i files..." % len(filenames))
        d = index_db(idx_filename, filenames, format)
        if len(d) == 0 and at_least_one_record(filenames, format):
            sys_exit("Index %s wrongly reports zero records" % idx_filename)
        print("%s - OK, %i records in %i files" %
              (idx_filename, len(d), len(filenames)))
    else:
        # One index per sequence file
        for filename in filenames:
            idx_filename = filename + ".idx"
            if options.reindex and os.path.isfile(idx_filename):
                print("%s - re-indexing..." % filename)
                os.remove(idx_filename)
            elif not os.path.isfile(idx_filename):
                print("%s - indexing..." % filename)
            d = index_db(idx_filename, filename, format)
            if len(d) == 0 and at_least_one_record([filename], format):
                sys_exit("Index %s wrongly reports zero records" %
                         idx_filename)
            print("%s - OK, %i records" % (idx_filename, len(d)))


if __name__ == "__main__":
    main()
