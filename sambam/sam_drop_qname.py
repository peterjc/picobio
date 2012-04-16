#!/usr/bin/env python
"""Python script to drop read name (QNAME) from SAM/BAM files.

This script is designed to be used as part of a Unix pipeline. It reads
SAM format data from stdin, and writes SAM format data to stdout.

The only change made to the SAM reads is in the QNAME field. For
single-fragment reads, QNAME is dropped (set to * for missing).
For multi-fragment reads (e.g. paired end reads), a QNAME is
required to cross reference the parts. Here short automatic names
are substituted instead.

The optional argument prefix is added to the start of any generated
read name (allowing you to avoid read name clashes).

Simple usage with SAM files:

$ ./sam_drop_names [prefix] < original.sam > dropped_names.sam

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_drop_names.py [prefix] | samtools view -S -b - > dropped_names.bam

If your SAM/BAM files lack @SQ headers, you may need to give
samtools the reference FASTA file as well.

Copyright Peter Cock 2012. All rights reserved. See:
https://github.com/peterjc/picobio
http://blastedbio.blogspot.co.uk/2012/03/bam-verus-cram-07.html
"""

import sys

if len(sys.argv) == 1:
    prefix = ""
elif len(sys.argv) == 2:
    prefix = sys.argv[1]
else:
    sys.stderr.write("Error, expect one optional parameter only (read name prefix)")
    sys.exit(1)

count = 0
mapping = dict()
#TODO - Automatically remove mapping entries once all parts of the read
#have been found? They would typically be near each other in the file...
#otherwise memory will be a problem with big paired end datasets.
for line in sys.stdin:
    if line[0]!="@":
        #Should be a read
        qname, flag, rest = line.split("\t", 2)
        if int(flag) & 0x1:
            #Multi-fragment read
            try:
                qname = prefix + str(mapping[qname])
            except KeyError:
                count += 1
                mapping[qname] = count
                qname = prefix + str(count)
        else:
            #Single fragment read
            qname = "*"
        line = "\t".join([qname, flag, rest])
    sys.stdout.write(line)
sys.stderr.write("Modified %i multi-fragment reads\n" % count)
