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
QUAL fields set to just *, so as part of the script it restores the SEQ
values.

Aim is to replace the combined effect of these three scripts:
* https://github.com/peterjc/picobio/blob/master/sambam/sam_restore_seq.py
* https://github.com/peterjc/picobio/blob/39fd9655e3641a22430c42b3d25fbcb4e6f4e099/blooming_reads/dedup_circular_sam.py
* https://github.com/peterjc/picobio/blob/master/blooming_reads/re_pair_circular_sam.py

TODO:

Review/use secondary alignment 0x100 FLAG when removing duplicates.

Review/set properly paired 0x2 FLAG after removing duplicates, and
recalculate TLEN. See also re_pair_circular_sam.py

Depending on feelings of SAM/BAM community, may want to split reads
spanning the origin into two sub-parts. Currently they spill over the
"end" of the sequence.

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys

from optparse import OptionParser


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)


VERSION = "0.0.1"

parser = OptionParser(usage="usage: %prog [options]\n\n" + usage,
                      version="%prog " + VERSION)

# References
parser.add_option("-l", "--lref", dest="linear_references",
                  type="string", metavar="FILE", action="append",
                  help="FASTA file of linear reference sequence(s). "
                       "Several files can be given if required.")
parser.add_option("-c", "--cref", dest="circular_references",
                  type="string", metavar="FILE", action="append",
                  help="FASTA file of circular reference sequence(s). "
                       "Several files can be given if required.")
# Reads
parser.add_option("-i", "--input", dest="input_reads",
                  type="string", metavar="FILE",
                  help="Input file of SAM format mapped reads to be processed (def. stdin)")
parser.add_option("-o", "--output", dest="output_reads",
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


def decode_cigar(cigar):
    """Returns a list of 2-tuples, integer count and operator char."""
    count = ""
    answer = []
    for letter in cigar:
        if letter.isdigit():
            count += letter  # string addition
        elif letter in "MIDNSHP=X":
            answer.append((int(count), letter))
            count = ""
        else:
            raise ValueError("Invalid character %s in CIGAR %s" %
                             (letter, cigar))
    return answer


assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [(
    14, 'S'), (15, 'M'), (1, 'P'), (1, 'D'), (3, 'P'), (54, 'M'), (1, 'D'), (34, 'M'), (5, 'S')]


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
    sys.stderr.write("Lengths of %i linear references loaded\n" %
                     len(ref_len_linear))
ref_len_circles = dict()
if options.circular_references:
    for f in options.circular_references:
        ref_len_circles.update(get_fasta_ids_and_lengths(f))
    sys.stderr.write("Lengths of %i circular references loaded\n" %
                     len(ref_len_circles))


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
                sys_exit("Bad SAM file, stay header lines?:\n%s%s" %
                         ("".join(batch), line))
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


def restore_seq(sam_lines):
    """Process a batch/list of SAM lines restoring missing SEQ.

    Designed for use on BWA MEM output where only the first read
    has the SEQ and QUAL recorded, additional alignments are just
    recorded with SEQ * and QUAL * instead.
    """
    global seq_mod
    last_name = None
    last_frag = None
    last_seq = None
    last_qual = None
    for line in sam_lines:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest = line.split(
            "\t", 11)
        if seq == "*":
            if qname == last_name and get_frag(flag) == last_frag:
                assert last_seq
                exp_len = cigar_seq_len(cigar)
                if exp_len < len(last_seq) and "H" in cigar:
                    # Ought to work if record it as soft trimming...
                    cigar = cigar.replace("H", "S")
                    exp_len = cigar_seq_len(cigar)
                assert exp_len == len(last_seq), \
                    "Cached SEQ %r length %i, but this read CIGAR expects length %i:\n%s" \
                    % (last_seq, len(last_seq), cigar_seq_len(cigar), line)
                seq = last_seq
                if qual == "*":
                    qual = last_qual
                seq_mod += 1
                line = "\t".join(
                    [qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest])
        elif "H" not in cigar:
            # Cache the SEQ
            last_name = qname
            last_frag = get_frag(flag)
            last_seq = seq
            last_qual = qual
        yield line


def undouble_circle_mappings(sam_lines):
    """Iterates over SAM lines updating POS for those mapped to circles."""
    global ref_len_circles
    for line in sam_lines:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, rest = line.split(
            "\t", 9)
        if rname in ref_len_circles:
            length = ref_len_circles[rname]
            int_pos = int(pos) - 1
            if int_pos != -1 and length <= int_pos:
                assert int_pos < length * 2, \
                    "Have POS %s yet length is %i or %i when doubled!\n%r" \
                    % (pos, length, length * 2, line)
                int_pos -= length
                pos = str(int_pos + 1)
                tlen = "0"  # old value invalidated
        if rnext in ref_len_circles:
            length = ref_len_circles[rnext]
            int_pnext = int(pnext) - 1
            if int_pnext != -1 and length <= int_pnext:
                assert int_pnext < length * 2, \
                    "Have PNEXT %s yet length is %i or %i when doubled!\n%r" \
                    % (pnext, length, length * 2, line)
                int_pnext -= length
                pnext = str(int_pnext + 1)
                tlen = "0"  # old value invalidated
        elif rnext == "=" and rname in ref_len_circles:
            length = ref_len_circles[rname]
            int_pnext = int(pnext) - 1
            if int_pnext != -1 and length <= int_pnext:
                assert int_pnext < length * 2, \
                    "Have PNEXT %s (%s) yet length is %i or %i when doubled!\n%r" \
                    % (pnext, rname, length, length * 2, line)
                int_pnext -= length
                pnext = str(int_pnext + 1)
                tlen = "0"  # old value invalidated
        yield "\t".join([qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, rest])


def dedup_batch(sam_lines):
    """De-duplicates and sorts a list of SAM lines for a single QNAME.

    Currently given a set of duplicates the first line is preserved.
    TODO - Look at the 0x100 FLAG for secondary alignment here?
    """
    global dup_reads_removed
    reads = {}
    batch_count = 0
    # sys.stderr.write("-"*80 + "\n")
    for line in sam_lines:
        batch_count += 1
        qname, flag, rname, pos, rest = line.split("\t", 4)
        int_flag = int(flag)
        frag = get_frag(int_flag)
        int_pos = int(pos) - 1
        rev_strand = bool(int_flag & 0x10)
        key = (qname, frag, rname, int_pos, rev_strand)
        if key in reads:
            dup_reads_removed += 1
            # sys.stderr.write("%s -Drop- %s" % (key, line))
        else:
            reads[key] = line
            # sys.stderr.write("%s -KEEP- %s" % (key, line))
    for key in sorted(reads):
        yield reads[key]
    # sys.stderr.write("%s %i --> %i\n" % (qname, batch_count, len(reads)))


# Open handles
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
seq_mod = 0
dup_reads_removed = 0
for batch in batch_by_qname(input_handle):
    # sys.stderr.write("%s\nBatch of %i lines:\n%s%s\n" % ("-" * 80, len(batch), "".join(batch), "-" * 80))
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
                if rname in ref_len_linear:
                    assert length == ref_len_linear[rname]
                    # print("Found @SQ line for linear reference %s" % rname)
                elif rname in ref_len_circles:
                    assert length == 2 * ref_len_circles[rname]
                    # Return the length to its correct value
                    # print("Fixing @SQ line for %s, length %i --> %i"
                    #       % (rname, length, ref_len_circles[rname]))
                    line = "@SQ\tSN:%s\tLN:%i\n" % (
                        rname, ref_len_circles[rname])
                elif rname is None:
                    sys_exit("Bad @SQ line:\n%s" % line)
                else:
                    sys_exit("This reference was not given!:\n%s" % line)
            output_handle.write(line)
        continue
    # Should be a batch of reads...
    for line in dedup_batch(undouble_circle_mappings(restore_seq(batch))):
        count += 1
        output_handle.write(line)

# Close handles
if options.input_reads:
    input_handle.close()
if options.output_reads:
    output_handle.close()

sys.stderr.write("Processed %i reads\n" % count)
sys.stderr.write("Restored missing SEQ for %i reads\n" % seq_mod)
sys.stderr.write(
    "Undoubling the circles removed %i duplicate alignments\n" % dup_reads_removed)
