#!/usr/bin/env python
usage = """Python script to remove SAM reads with long CIGAR strings.

The BAM format (currently) uses an unsigned 16bit integer for the
number of CIGAR operations in a read, and therefore BAM files can
only hold reads with up to 65535 CIGAR operators. SAM does not have
this limit, but the samtools implementation (reasonably) also has
the same 16bit limit. See also:
https://github.com/samtools/samtools/pull/39

This script is designed to be used as part of a Unix pipeline. It
takes no command line arguments. It reads SAM format data from stdin,
and writes SAM format data to stdout.

The only change made to the SAM reads is to drop records with over
65535 CIGAR operators. These are logged to stderr.

$ ./sam_drop_long_cigar.py < original.sam > no_long_cigar.sam

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_drop_long_cigar.py | samtools view -S -b - > no_long_cigar.bam

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys

if len(sys.argv) != 1:
    sys.stderr.write("ERROR: Bad arguments.\n\n")
    sys.stderr.write("Expects SAM on stdin, and writes SAM to stdout.\n")
    sys.exit(1)

#def decode_cigar(cigar):
#    """Returns a list of 2-tuples, integer count and operator char."""
#    count = ""
#    answer = []
#    for letter in cigar:
#        if letter.isdigit():
#            count += letter #string addition
#        elif letter in "MIDNSHP=X":
#            answer.append((int(count), letter))
#            count = ""
#        else:
#            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
#    return answer
#
#assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [(14,'S'),(15,'M'),(1,'P'),(1,'D'),(3,'P'),(54,'M'),(1,'D'),(34,'M'),(5,'S')]

def cigar_length(cigar):
    """Returns number of cigar operators (integer)."""
    answer = 0
    for letter in cigar:
        if letter.isdigit():
            pass
        elif letter in "MIDNSHP=X":
            answer += 1
        else:
            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
    return answer

count = 0
longs = 0
for line in sys.stdin:
    if line[0]!="@":
        #Should be a read
        count += 1
        qname, flag, rname, pos, mapq, cigar, rest = line.split("\t", 6)
        if cigar != "*":
            len_cigar = cigar_length(cigar)
            if len_cigar > 65535:
                longs += 1
                sys.stderr.write("Dropping read %s with %i CIGAR operators\n" % (qname, len_cigar))
                continue
    sys.stdout.write(line)
sys.stderr.write("Dropped %i out of %i reads\n" % (longs, count))
