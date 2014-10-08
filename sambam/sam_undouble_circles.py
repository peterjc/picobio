#!/usr/bin/env python
usage = """Post-process SAM mappings onto doubled circular references.

Most read mappers do not support cicular reference sequences. One
workarround is to make a doubled FASTA file and map onto that, then
post process to ensure POS is in the range 1 to circle length, and
remove any artifact duplicated mappings.

This script is designed to be used as part of a Unix pipeline. It reads
SAM format data from stdin, and writes SAM format data to stdout.

This script is designed to be used on the BWA-MEM output where it seems
with the -a option additional alignments are reported with the SEQ and
QUAL fields set to just *, so as part of the script is restores the SEQ
values.

Aim is to replace the combined effect of these three scripts:
* https://github.com/peterjc/picobio/blob/master/sambam/sam_restore_seq.py
* https://github.com/peterjc/picobio/blob/master/blooming_reads/dedup_circular_sam.py
* https://github.com/peterjc/picobio/blob/master/blooming_reads/re_pair_circular_sam.py

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys
import os
from optparse import OptionParser

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

VERSION = "0.0.0"

parser = OptionParser(usage="usage: %prog [options]",
                      version="%prog "+VERSION)

#References
parser.add_option("-l", "--lref", dest="linear_references",
                  type="string", metavar="FILE", action="append",
                  help="FASTA file of linear reference sequence(s). "
                       "Several files can be given if required.")
parser.add_option("-c", "--cref", dest="circular_references",
                  type="string", metavar="FILE", action="append",
                  help="FASTA file of circular reference sequence(s). "
                       "Several files can be given if required.")
#Reads
parser.add_option("-i", "--input", dest="input_reads",
                  type="string", metavar="FILE",
                  help="Input file of SAM format mapped reads to be processed (def. stdin)")
parser.add_option("-o","--output", dest="output_reads",
                  type="string", metavar="FILE",
                  help="Output file for processed SAM format mapping (def. stdout)")

(options, args) = parser.parse_args()

if len(sys.argv)==1:
    parser.print_help()
    sys.exit(1)

if (not options.linear_references) and (not options.circular_references):
    parser.error("You must supply some linear and/or circular references")
if args:
    parser.error("No arguments expected")


def decode_cigar(cigar):
    """Returns a list of 2-tuples, integer count and operator char."""
    count = ""
    answer = []
    for letter in cigar:
        if letter.isdigit():
            count += letter #string addition
        elif letter in "MIDNSHP=X":
            answer.append((int(count), letter))
            count = ""
        else:
            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
    return answer

assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [(14,'S'),(15,'M'),(1,'P'),(1,'D'),(3,'P'),(54,'M'),(1,'D'),(34,'M'),(5,'S')]

def cigar_seq_len(cigar_str):
    slen = 0
    for count, operator in decode_cigar(cigar_str):
        if operator in "MIS=X":
            slen += count
    return slen
assert cigar_seq_len("1I58M1I34M1I2M1D12M") == 109

def get_frag(flag):
    f = int(flag)
    if f & 0x1:
        #multi-part
        first = f & 0x40
        last = f & 0x80
        if first and last:
            return None # Part of a mult-fragment read, not pair?
        elif first:
            return 1
        elif last:
            return 2
        else:
            return None # Unknown
    else:
        return None



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


#Load the reference sequence lengths
ref_len_linear = dict()
if options.linear_references:
    for f in options.linear_references:
        ref_len_linear.update(get_fasta_ids_and_lengths(f))
    sys.stderr.write("Lengths of %i linear references loaded\n" % len(ref_len_linear))
ref_len_circles = dict()
if options.circular_references:
    for f in options.circular_references:
        ref_len_circles.update(get_fasta_ids_and_lengths(f))
    sys.stderr.write("Lengths of %i circular references loaded\n" % len(ref_len_circles))


#Open handles
if options.input_reads:
    input_handle = open(options.input_reads)
else:
    input_handle = sys.stdin
if options.output_reads:
    output_handle = open(options.output_reads, "w")
else:
    output_handle = sys.stdout

last_seq = None
last_qual = None
last_name = None
last_frag = 0
count = 0
mod = 0
for line in input_handle:
    if line[0]!="@":
        #Should be a read
        count += 1
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest = line.split("\t", 11)
        if seq == "*":
            if qname == last_name and get_frag(flag) == last_frag:
                assert last_seq
                exp_len = cigar_seq_len(cigar)
                if exp_len < len(last_seq) and "H" in cigar:
                    # Ought to work if record it as sort trimming...
                    cigar = cigar.replace("H", "S")
                    exp_len = cigar_seq_len(cigar)
                assert exp_len == len(last_seq), \
                    "Cached SEQ %r length %i, but this read CIGAR expects length %i:\n%s" \
                    % (last_seq, len(last_seq), cigar_seq_len(cigar), line)
                seq = last_seq
                if qual == "*":
                    qual = last_qual
                mod += 1
                line = "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest])
        elif "H" not in cigar:
            #Cache the SEQ
            last_name = qname
            last_frag = get_frag(flag)
            last_seq = seq
            last_qual = qual
    output_handle.write(line)

#Close handles
if options.input_reads:
    input_handle.close()
if options.output_reads:
    output_handle.close()

sys.stderr.write("Modified %i out of %i reads\n" % (mod, count))
