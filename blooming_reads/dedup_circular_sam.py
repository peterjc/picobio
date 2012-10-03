#!/usr/bin/env python
"""Python script to remap sequencing reads to a circular genomes.

The idea is that you take your circular genomes, and double them,
and use a traditional mapper on that (in single end mode!).
e.g. I'm mrfast at the moment for this,

What this script does first is remove duplicates due to mapping
in the second half (or should shift them to the first half if
just getting one position per read).

TODO:

Step two would be to split reads mapped over the origin into
two parts (to comply with the SAM/BAM specification). Without
doing this it seems Tablet v1.12.09.03 will show spill over in
SAM, but crops the view for BAM.

Step three would be to review paired end information, and take
into consideration the circular nature when deciding if a given
pair mapping is sensible, and if you have multiple mappings
which is the best,

One of my aims here is to explore how best to output this kind of
mapping in SAM/BAM format (currently not defined), which will
probably mean splitting a read mapping over the origin into two
fragments (i.e. two lines in SAM). This means that with paired
end data, you might get two, three or even four lines in SAM
(rather than the normal two lines, one for each half of the pair).
"""

import sys
import os
import tempfile
import time
from optparse import OptionParser

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

VERSION = "0.0.0"


def go(input, output, paired, linear_refs, circular_refs):
    ref_len_linear = dict()
    if linear_refs:
        for f in linear_refs:
            ref_len_linear.update(get_fasta_ids_and_lengths(f))
    ref_len_circles = dict()
    if circular_refs:
        for f in circular_refs:
            ref_len_circles.update(get_fasta_ids_and_lengths(f))
    #print ref_len_circles

    if input is None:
        input_handle = sys.stdin
    elif isinstance(input, basestring):
        input_handle = open(input)
    else:
        input_handle = input

    if output is None:
        output_handle = sys.stdout
    elif isinstance(output, basestring):
        output_handle = open(output)
    else:
        output_handle = output

    line = input_handle.readline()
    while line[0] == "@":
        #SAM header
        if line[0:4] == "@SQ\t":
            parts = line[4:].strip().split("\t")
            rname = None
            length = None
            for p in parts:
                if p.startswith("SN:"):
                    rname = p[3:]
                if p.startswith("LN:"):
                    length = int(p[3:])
            if rname in ref_len_linear:
                assert length == ref_len_linear[rname]
                #print "Found @SQ line for linear reference %s" % rname
            elif rname in ref_len_circles:
                assert length == 2 * ref_len_circles[rname]
                #Return the length to its correct value
                #print "Fixing @SQ line for %s, length %i --> %i" % (rname, length, ref_len_circles[rname])
                line = "@SQ\tSN:%s\tLN:%i\n" % (rname, ref_len_circles[rname])
            elif rname is None:
                sys_exit("Bad @SQ line:\n%s" % line)
            else:
                sys_exit("This reference was not given!:\n%s" % line)
        output_handle.write(line)
        line = input_handle.readline()

    prev_line = None
    while line:
        #SAM read
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, rest = line.split("\t", 10)
        if rname in ref_len_circles and pos != "0":
            length = ref_len_circles[rname]
            pos = int(pos) - 1
            if length <= pos:
                assert pos < length*2, "Have POS %i yet length is %i or %i when doubled!\n%r" % (pos, length, length*2, line)
                pos -= length #Modulo circle length
                pos = str(pos+1)
                line = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, rest])
                if line == prev_line:
                    #print "Dropping duplicate"
                    line = ""
                else:
                    #print "Fixed: %r" % line
                    pass
        #Next line...
        output_handle.write(line)
        prev_line = line
        line = input_handle.readline()


    if isinstance(input, basestring):
        input_handle.close()
    if isinstance(output, basestring):
        output_handle.close()


def get_fasta_ids_and_lengths(fasta_filename):
    h = open(fasta_filename)
    name = None
    length = 0
    for line in h:
        if line[0] == ">":
            if name is not None:
                yield name, length
            name = line[1:].split(None, 1)[0]
            length = 0
        else:
            length += len(line.strip())
    if name is not None:
        yield name, length
    h.close()

def main():
    parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog "+VERSION)
    #References
    parser.add_option("-l", "--lref", dest="linear_references",
                      type="string", metavar="FILE", action="append",
                      help="""FASTA file of linear reference sequence(s)
                           Several files can be given if required.""")
    parser.add_option("-c", "--cref", dest="circular_references",
                      type="string", metavar="FILE", action="append",
                      help="""FASTA file of circular reference sequence(s)
                           Several files can be given if required.""")
    #Reads
    #TODO - Make paired mode or single mode the default?
    parser.add_option("-i", "--input", dest="input_reads",
                      type="string", metavar="FILE",
                      help="Input file of SAM format mapped reads to be processed (def. stdin)")
    parser.add_option("-o","--output", dest="output_reads",
                      type="string", metavar="FILE",
                      help="Output file for processed SAM format mapping (def. stdout)")
    
    (options, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if (not options.linear_references) and (not options.circular_references):
        parser.error("You must supply some linear and/or circular references")

    if args:
        parser.error("No arguments expected")

    paired = True
    go(options.input_reads, options.output_reads, paired,
       options.linear_references, options.circular_references)

if __name__ == "__main__":
    main()
