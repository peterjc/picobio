#!/usr/bin/env python
usage = """Python script to remove paired information in SAM reads.

The intended usage is where you wish to treat "orphaned" paired
reads as single reads, meaning removing any /1 or /2 suffix in
the FASTQ file and likewise clearing the paired bits in the SAM
FLAG.

This script is designed to be used as part of a Unix pipeline. It
takes no command line arguments. It reads SAM format data from stdin,
and writes SAM format data to stdout.

The only change made to the FLAG field, clearing the following bits:
* 0x1 template having multiple segments in sequencing
* 0x8 next segment in the template unmapped
* 0x20 next segment mapped to reverse strand
* 0x40 the first segment in the template
* 0x80 the last segment in the template

Example:

$ ./sam_depair.py < original.sam > as_singles.sam

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_depair.py | samtools view -S -b - > as_singles.bam

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys

if len(sys.argv) != 1:
    sys.stderr.write("ERROR: Bad arguments.\n\n")
    sys.stderr.write("Expects SAM on stdin, and writes SAM to stdout.\n")
    sys.exit(1)

count = 0
tweaked = 0
mask = 0x1 | 0x8 | 0x20 | 0x40 | 0x80
flip_mask = ~mask
for line in sys.stdin:
    if line[0] != "@":
        # Should be a read
        count += 1
        qname, flag, rest = line.split("\t", 2)
        flag = int(flag)
        if flag & mask:
            # Want to clear those bits...
            flag = flag & flip_mask
            tweaked += 1
            line = "%s\t%i\t%s" % (qname, flag, rest)
    sys.stdout.write(line)
sys.stderr.write("Tweaked %i out of %i reads\n" % (tweaked, count))
