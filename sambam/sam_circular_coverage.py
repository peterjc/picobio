#!/usr/bin/env python
usage = """Calculate SAM read coverage allowing for circular genomes.

Assumes your SAM file where reads mapping over the origin of a
circular reference spill over (aligned end exceeds length), or
are directly from mapping against a doubled-reference - and thus
calculates coverage at each position modulo the circle length.

Currently calculates coverage split into these five categories:
- Single reads all mapped to same reference
- Single reads mapped to multiple references
- Paired reads both mapped to same reference
- Paired reads mapped to multiple references
- Paired reads with partner unmapped

The output is a FASTA QUAL like plain text file, ">identifier"
followed by five lines of space separated coverage scores (one
value for each base in the associated reference sequence). If
all the entries for a category are zero "None" is record instead.

TODO: Switch to JSON output?
"""

import sys

from optparse import OptionParser
from builtins import range  # for Python 2

import numpy as np


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)


VERSION = "0.0.3"

parser = OptionParser(
    usage="usage: %prog [options]\n\n" + usage, version="%prog " + VERSION
)
# References
parser.add_option(
    "-l",
    "--lref",
    dest="linear_references",
    type="string",
    metavar="FILE",
    action="append",
    help="FASTA file of linear reference sequence(s). "
    "Several files can be given if required.",
)
parser.add_option(
    "-c",
    "--cref",
    dest="circular_references",
    type="string",
    metavar="FILE",
    action="append",
    help="FASTA file of circular reference sequence(s). "
    "Several files can be given if required.",
)
# Reads
parser.add_option(
    "-i",
    "--input",
    dest="input_reads",
    type="string",
    metavar="FILE",
    help="Input file of SAM format mapped reads to be processed (def. stdin)",
)
parser.add_option(
    "-o",
    "--output",
    dest="coverage_file",
    type="string",
    metavar="FILE",
    help="Output file for coverage report (def. stdout)",
)

(options, args) = parser.parse_args()

if len(sys.argv) == 1:
    parser.print_help()
    sys.exit(1)

if (not options.linear_references) and (not options.circular_references):
    parser.error("You must supply some linear and/or circular references")
if args:
    parser.error("No arguments expected")

solo0 = solo1 = solo2 = solo12 = 0


def go(input_handle, output_handle, linear_refs, circular_refs):
    sam_len_references = dict()

    # Should be a batch of reads...
    global solo0, solo1, solo2, solo12
    solo0 = solo1 = solo2 = solo12 = 0

    global coverage
    coverage = dict()

    for lengths in [ref_len_linear, ref_len_circles]:
        for ref, length in lengths.iteritems():
            coverage[ref] = np.zeros((5, length), np.float)

    for batch in batch_by_qname(input_handle):
        if not batch:
            continue
        if batch[0][0] == "@":
            # SAM header
            for line in batch:
                assert line[0] == "@"
                if line[0:4] == "@SQ\t":
                    parts = line[4:].strip().split("\t")
                    rname = None
                    length = None
                    for p in parts:
                        if p.startswith("SN:"):
                            rname = p[3:]
                        if p.startswith("LN:"):
                            length = int(p[3:])
                    sam_len_references[rname] = length
                    if rname in ref_len_linear:
                        assert length == ref_len_linear[rname]
                        # print("Found @SQ line for linear reference %s"
                        #       % rname)
                    elif rname in ref_len_circles:
                        if length == 2 * ref_len_circles[rname]:
                            # We will use the length from the FASTA file
                            sys.stderr.write(
                                "WARNING: @SQ line for %s gives length %i, double %i in FASTA\n"
                                % (rname, length, ref_len_circles[rname])
                            )
                        else:
                            assert length == ref_len_circles[rname]
                    elif rname is None or length is None:
                        sys_exit("Bad @SQ line:\n%s" % line)
                    else:
                        sys_exit("This reference was not given!:\n%s" % line)
            # End of header
            continue

        # Split up all the QNAME reads by fragment
        r0 = set()
        r1 = set()
        r2 = set()
        rnames = set()
        reads = {None: r0, 1: r1, 2: r2}
        qname = None
        for line in batch:
            # SAM read
            qname, flag, rname, pos, mapq, cigar, rest = line.split("\t", 6)
            flag = int(flag)
            if flag & 0x4:
                # Unmapped, ignore
                continue
            frag = get_frag(flag)
            reads[frag].add((rname, pos, cigar))
            rnames.add(rname)  # to see if all map to same ref
        # Now get weighted coverage
        if r0:
            # This QNAME is a single read (perhaps multiply mapped)
            assert not r1 and not r2, (
                "Inconsistent FLAG for %s, is it paired or unpaired?" % qname
            )
            solo0 += 1
            if len(rnames) == 1:
                field = 0
            else:
                field = 1
            weight = 1.0 / len(r0)
            for rname, pos, cigar in r0:
                count_coverage(coverage, field, weight, rname, pos, cigar)
        elif r1 or r2:
            # This QNAME is a paired read
            if r1 and r2:
                solo12 += 1
                if len(rnames) == 1:
                    field = 2
                else:
                    field = 3
            elif r1:
                solo1 += 1
                field = 4
            elif r2:
                solo2 += 1
                field = 4
            else:
                assert False, "Error sorting %s by fragment" % qname
            if r1:
                weight = 1.0 / len(r1)
                for rname, pos, cigar in r1:
                    count_coverage(coverage, field, weight, rname, pos, cigar)
            if r2:
                weight = 1.0 / len(r2)
                for rname, pos, cigar in r2:
                    count_coverage(coverage, field, weight, rname, pos, cigar)

    for lengths in [ref_len_linear, ref_len_circles]:
        for ref, length in lengths.iteritems():
            output_handle.write(">%s length %i\n" % (ref, length))
            m = 0
            for row in coverage[ref]:
                assert len(row) == length
                row_max = max(row)
                if row_max:
                    output_handle.write("\t".join("%.1f" % v for v in row) + "\n")
                else:
                    output_handle.write("None\n")
                m = max(m, row_max)
            sys.stderr.write("%s length %i max depth %i\n" % (ref, length, m))
    sys.stderr.write(
        "%i singletons; %i where only /1, %i where only /2, %i where both mapped\n"
        % (solo0, solo1, solo2, solo12)
    )


def cigar_tuples(cigar_str):
    """CIGAR string parsed into a list of tuples (operator code, count).

    e.g. cigar string of 36M2I3M becomes [('M', 36), ('I', 2), ('M', 3)]

    Any empty CIGAR string (represented as * in SAM) is given as None.
    """
    if cigar_str == "*":
        return None
    answer = []
    count = ""
    for letter in cigar_str:
        if letter.isdigit():
            count += letter  # string addition
        else:
            if letter not in "MIDNSHP=X":
                raise ValueError(
                    "Invalid character %s in CIGAR %s" % (letter, cigar_str)
                )
            answer.append((letter, int(count)))
            count = ""
    return answer


def cigar_alen(cigar_str):
    alen = 0
    for operator, count in cigar_tuples(cigar_str):
        if operator in "MDN=X":
            alen += count
    return alen


def count_coverage(coverage, field, weight, rname, pos, cigar):
    pos = int(pos) - 1
    values = coverage[rname]
    length = values.shape[1]
    for i in range(pos, pos + cigar_alen(cigar)):
        values[field, i % length] += weight


def get_frag(flag):
    f = int(flag)
    if f & 0x1:
        # multi-part
        first = f & 0x40
        last = f & 0x80
        if first and last:
            return None  # Part of a mult-fragment read, not pair?
        elif first:
            return 1
        elif last:
            return 2
        else:
            return None  # Unknown
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


# Load the reference sequence lengths
ref_len_linear = dict()
if options.linear_references:
    for f in options.linear_references:
        ref_len_linear.update(get_fasta_ids_and_lengths(f))
    sys.stderr.write("Lengths of %i linear references loaded\n" % len(ref_len_linear))
ref_len_circles = dict()
if options.circular_references:
    for f in options.circular_references:
        ref_len_circles.update(get_fasta_ids_and_lengths(f))
    sys.stderr.write(
        "Lengths of %i circular references loaded\n" % len(ref_len_circles)
    )


def batch_by_qname(input_handle):
    """Yields lists of SAM lines, batching by read name.

    If there is a SAM header, that is returned first.

    There after you get all the SAM lines for read name one, then
    for read name two, etc. This assumes the SAM file is sorted
    or at least grouped by read name (typical of alignment output).
    """
    batch = []
    batch_qname = None
    for line in input_handle:
        if line[0] == "@":
            # SAM Header
            if batch_qname or (batch and batch[-1][0] != "@"):
                sys_exit(
                    "Bad SAM file, stay header lines?:\n%s%s" % ("".join(batch), line)
                )
            batch.append(line)
        else:
            # SAM read
            qname, rest = line.split("\t", 1)
            if batch_qname == qname:
                batch.append(line)
            else:
                yield batch
                batch = [line]
                batch_qname = qname
    # End of file
    if batch:
        yield batch


# Open handles
if options.input_reads:
    input_handle = open(options.input_reads)
else:
    input_handle = sys.stdin
if options.coverage_file:
    output_handle = open(options.coverage_file, "w")
else:
    output_handle = sys.stdout


go(input_handle, output_handle, ref_len_circles, ref_len_circles)


# Close handles
if options.input_reads:
    input_handle.close()
if options.coverage_file:
    output_handle.close()
