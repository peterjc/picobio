#!/usr/bin/env python
"""Simple FASTQ interlacer.

Checks read identifiers agree, or end with /1 and /2 respectively.
"""
import sys
import gzip
from Bio.SeqIO.QualityIO import FastqGeneralIterator

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

if len(sys.argv) != 3:
    sys_exit("Requires two arguments, a pair of FASTQ filenames")
fastq1 = sys.argv[1]
fastq2 = sys.argv[2]

sys.stderr.write("Interlacing %s and %s\n" % (fastq1, fastq2))
if fastq1.endswith(".gz"):
    sys.stderr.write("Decompressing %s\n" % fastq1)
    handle1 = gzip.open(fastq1)
else:
    handle1 = open(fastq1)
if fastq2.endswith(".gz"):
    sys.stderr.write("Decompressing %s\n" % fastq2)
    handle2 = gzip.open(fastq2)
else:
    handle2 = open(fastq2)
sys.stderr.write("Interlacing paired FASTQ files to stdout...\n")
out_handle = sys.stdout

iter1 = FastqGeneralIterator(handle1)
iter2 = FastqGeneralIterator(handle2)

for title1, seq1, qual1 in iter1:
    try:
        title2, seq2, qual2 = iter2.next()
    except StopIteration:
        sys_exit("More records in %s than %s, e.g. %s" % (fastq1, fastq2, title1))
    id1, descr1 = title1.split(None, 1)
    id2, descr2 = title2.split(None, 1)
    if id1 == id2:
        #Add the /1 and /2, preserve any description after the ID
        if descr1:
            descr1 = " " + descr1
        if descr2:
            descr2 = " " + descr2
        out_handle.write("@%s/1%s\n%s\n+\n%s\n@%s/2%s\n%s\n+\n%s\n" \
                         % (id1, descr1, seq1, qual1, id2, descr2, seq2, qual2))
    elif id1.endswith("/1") and id2.endswith("/2") and id1[:-2]==id2[:-2]:
        out_handle.write("@%s\n%s\n+\n%s\n@%s\n%s\n+\n%s\n" \
                         % (title1, seq1, qual1, title2, seq2, qual2))
    else:
        sys_exit("Mismatched records %r vs %r" % (title1, title2))

#Check at end of file two
try:
    title2,seq2, qual2 = iter2.next()
    sys_exit("More records in %s than %s, e.g. %s" % (fastq2, fastq1, title2))
except StopIteration:
    pass

handle1.close()
handle2.close()
sys.stderr.write("Interlacing paired FASTQ files done.\n")
