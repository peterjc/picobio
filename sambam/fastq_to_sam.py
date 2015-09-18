#!/usr/bin/env python
"""Python script to turn FASTQ into unaliged SAM/BAM files.

This script is designed to be used as part of a Unix pipeline. It
works with Python 2 and Python 3, e.g.

    $ python fastq_to_sam.py R1.fastq R2.fastq > unmapped.sam
    Done, 532 pairs

Or:

    $ python3 fastq_to_sam.py R1.fastq R2.fastq > unmapped.sam
    Done, 532 pairs

As long as the Python script is marked as executable you can do:

    $ ./fastq_to_sam.py R1.fastq R2.fastq > unmapped.sam
    Done, 532 pairs

Simple usage with BAM files with conversion to/from SAM via samtools:

    $ ./fastq_to_sam.py R1.fastq R2.fastq | samtools view -S -b - > unmapped.bam
    [samopen] no @SQ lines in the header.
    Done, 532 pairs

Note that no @SQ lines are expected in SAM/BAM files with only unaligned reads.

WARNING: This assumes your FASTQ files use the Sanger quality encoding.

TODO:

 - Test cases
 - Galaxy wrapper?
 - Proper command line API
 - Support for gzipped FASTQ (detected via filename?)
 - Support for interlaced FASTQ
 - Support for setting read groups
 - Support for multiple FASTQ input pairs (and read groups)

Copyright Peter Cock 2015. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys

if "-v" in sys.argv or "--version" in sys.argv:
    print("This is fastq_to_sam.py version 0.0.1")
    sys.exit(0)

# TODO - proper API, allow interleaved FASTQ, read group, etc
if len(sys.argv) != 3:
    sys.stderr.write("Expects two arguments, a pair of FASTQ filenames\n")
    sys.exit(1)

try:
    from Bio._py3k import zip
    from Bio.SeqIO.QualityIO import FastqGeneralIterator
except ImportError:
    sys.exit("ERROR: This requires Biopython.\n")
    sys.exit(1)

fastq1 = FastqGeneralIterator(open(sys.argv[1]))
fastq2 = FastqGeneralIterator(open(sys.argv[2]))

# Paired, unmapped, mate unmapped, either first or second in pair:
flag1 = "77"
flag2 = "141"
rname = "*"
pos = "0"
mapq = "0"
cigar = "*"
rnext = "*"
pnext = "0"
tlen = "0"

pairs = 0
for (t1, s1, q1), (t2, s2, q2) in zip(fastq1, fastq2):
    id1 = t1.split(None, 1)[0]
    id2 = t2.split(None, 1)[0]
    if id1 == id2:
        # Good, should we check the description follow Illumina naming?
        qname = id1
    else:
        assert id1.endswith("/1"), t1
        assert id2.endswith("/2"), t2
        qname = id1[:-2]
    
    print("\t".join([qname, flag1, rname, pos, mapq, cigar, rnext, pnext, tlen, s1, q1]))
    print("\t".join([qname, flag1, rname, pos, mapq, cigar, rnext, pnext, tlen, s2, q2]))
    pairs += 1
sys.stderr.write("Done, %i pairs\n" % pairs)
