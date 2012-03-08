#!/usr/bin/env python
"""Python script to remove tags from SAM/BAM files.

This script is designed to be used as part of a Unix pipeline. It
takes as optional command line arguments a white list of tags to
preserve. It reads SAM format data from stdin, and writes SAM format
data to stdout.

The only change made to the SAM reads is in the SEQ field of mapped
reads. Any bases matching the reference are replaced with equals
signs. This makes the files much easier to compress, which can be
demonstrated by comparing a gzipped version of the SAM files with
and without the equals, or their BAM equivalents.

Simple usage with SAM files, keeping only read-group tags:

$ ./sam_strip_tags.py RG < original.sam > equals.sam

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_strip_tags.py RG | samtools view -S -b - > equals.bam

If your SAM/BAM files lack @SQ headers, you may need to give
samtools the reference FASTA file as well.
"""

import sys

white_list = set(x.strip() for x in sys.argv[1:])

for line in sys.stdin:
    if line[0]!="@":
        #Should be a read
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags = line.rstrip().split("\t", 11)
        tags = "\t".join(t for t in tags.split("\t") if t[:2] in white_list)
        line = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags]) + "\n"
    sys.stdout.write(line)
