#!/usr/bin/env python
"""Python script to check if BGZF (e.g. BAM) files have EOF marker.

BAM files are compressed using BGZF, Blocked GNU Zip Format, which
is a variant of GZIP. Modern BAM files include a special empty
block at the end of the file (EOF) as a marker to help spot when
a dataset has been truncated. This is just a 28 byte BGZF block,
which when decompressed is empty.

Some early tools output valid BAM files without this optional
(but recommended) EOF marker.

Usage with one or more BAM or BGZF files:

$ ./bgzf_check_eof.py example1.bam example2.bam ... exampleN.bam

The filenames are checked in the order given, if any are invalid
the tool exits with a non-zero error level and a message to stderr.
If all the files are valid, it returns with a zero error level.

Return codes:
* 0 - No errors found
* 1 - Invalid arguments
* 2 - File not found
* 3 - File is zero bytes (and thus not valid BGZF or BAM)
* 4 - File missing BGZF header
* 5 - File looks like BGZF, but missing BGZF EOF marker

See also: http://samtools.sourceforge.net/

v0.0.0 - Original script
v0.0.1 - Dropped internal function sys_exit
"""

import os
import sys


def check_bam(filename):
    header = "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00\x42\x43\x02\x00"
    eof = (
        "\x1f\x8b\x08\x04\x00\x00\x00\x00\x00\xff\x06\x00BC"
        "\x02\x00\x1b\x00\x03\x00\x00\x00\x00\x00\x00\x00\x00\x00"
    )
    if not os.path.isfile(filename):
        sys.stderr.write("Missing file %s\n" % filename)
        sys.exit(2)
    size = os.path.getsize(filename)
    if not size:
        sys.stderr.write("Empty file (zero bytes) %s\n" % filename)
        sys.exit(3)
    h = open(filename, "rb")
    # Check it looks like a BGZF file
    # (could still be GZIP'd, in which case the extra block is harmless)
    data = h.read(len(header))
    if data != header:
        sys.stderr.write("File %s is not a BGZF/BAM file\n" % filename)
        sys.exit(4)
    # Check if it has the EOF already
    h.seek(size - 28)
    data = h.read(28)
    h.close()
    if data == eof:
        sys.stderr.write("Good, BGZF EOF already present in %s\n" % filename)
    else:
        sys.stderr.write("Missing EOF marker in BGZF/BAM file %s\n" % filename)
        sys.exit(5)


if len(sys.argv) == 1:
    sys.stderr.write(
        "Takes one or more BGZF/BAM filenames as arguments (edits in place)"
    )
    sys.exit(1)
for bam_filename in sys.argv[1:]:
    check_bam(bam_filename)
