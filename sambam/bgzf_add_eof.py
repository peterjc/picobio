#!/usr/bin/env python
"""Python script to add missing EOF marker to BAM or BGZF files.

BAM files are compressed using BGZF, Blocked GNU Zip Format, which
is a variant of GZIP. Modern BAM files include a special empty
block at the end of the file (EOF) as a marker to help spot when
a dataset has been truncated. This is just a 28 byte BGZF block,
which when decompressed is empty.

Some early tools output valid BAM files without this optional
(but recommended) EOF marker.

This script will add the EOF marker is not already present.

WARNING: If your BAM or BGZF file is truely truncated, this will
not magically fix it. It may hide or obscure the true problem.

WARNING: To avoid excessive data writing, this script modifies
the BAM or BGZF file in situ!

Usage with one or more BAM or BGZF files:

$ ./bam_add_eof.py example1.bam example2.bam ... exampleN.bam

See also: http://samtools.sourceforge.net/

v0.0.0 - Original script
v0.0.1 - Use append mode to add EOF block

"""

import os
import sys

def sys_exit(msg, return_code=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(return_code)

def fix_bam(filename):
    header = "\x1f\x8b\x08\x04\x00\x00\x00\x00" + \
             "\x00\xff\x06\x00\x42\x43\x02\x00"
    eof = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC" + \
          "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    if not os.path.isfile(filename):
        sys_exit("Missing file %s" % filename)
    size = os.path.getsize(filename)
    h = open(filename, "rb") #read only for now
    #Check it looks like a BGZF file
    #(could still be GZIP'd, in which case the extra block is harmless)
    data = h.read(len(header))
    if data != header:
        sys_exit("File %s is not a BAM file" % filename)
    #Check if it has the EOF already
    h.seek(size - 28)
    data = h.read(28)
    h.close()
    if data == eof:
        sys.stderr.write("EOF already present in %s\n" % filename)
    else:
        sys.stderr.write("Adding EOF block to %s\n" % filename)
        h = open(filename, "ab")
        h.write(eof)
        h.close()

if len(sys.argv) == 1:
    sys_exit("Takes one or more BGZF/BAM filenames as arguments (edits in place)")
for bam_filename in sys.argv[1:]:
    fix_bam(bam_filename)
