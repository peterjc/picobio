#!/usr/bin/env python
usage = """Python script to demonstrate reference based compression in SAM/BAM.

This script is designed to be used as part of a Unix pipeline. It
takes a single command line argument of a FASTA reference genome
filename. It reads SAM format data from stdin, and writes SAM format
data to stdout.

The only change made to the SAM reads is in the SEQ field of mapped
reads. Any bases matching the reference are replaced with equals
signs. This makes the files much easier to compress, which can be
demonstrated by comparing a gzipped version of the SAM files with
and without the equals, or their BAM equivalents.

Simple usage with SAM files, optional argument [mode] can be "add"
(default, insert equals signs) or "remove" (remove equals signs),
or "full" (set SEQ and QUAL to "*" for perfect matches and uses
the CIGAR equal operator, otherwise like the default "add" mode,
just add equals signs to SEQ for matching bases):

$ ./sam_seq_equals reference.fasta [mode] < original.sam > equals.sam

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_seq_equals reference.fasta [mode] | samtools view -S -b - > equals.bam

If your SAM/BAM files lack @SQ headers, you may need to give
samtools the reference FASTA file as well.

Copyright Peter Cock 2012. All rights reserved. See:
https://github.com/peterjc/picobio
http://blastedbio.blogspot.co.uk/2012/02/reference-based-sambam-compression.html
"""

import sys

if len(sys.argv) == 2:
    reference_filename = sys.argv[1]
    add_equals = True
    drop_seq = False
elif len(sys.argv) == 3:
    reference_filename = sys.argv[1]
    if sys.argv[2].lower() == "add":
        add_equals = True
        drop_seq = False
    elif sys.argv[2].lower() == "remove":
        add_equals = False
        drop_seq = False
    elif sys.argv[2].lower() == "full":
        add_equals = True
        drop_seq = True
    else:
        sys.stderr.write(
            "ERROR: Second (optional) argument must be 'add' (default) or 'remove' (no quotes)\n\n")
        sys.stderr.write(usage)
        sys.exit(1)
else:
    sys.stderr.write("ERROR: Bad arguments.\n\n")
    sys.stderr.write(
        "Expects SAM on stdin, and writes SAM to stdout. Requires the FASTA\n")
    sys.stderr.write(
        "reference filename as the first argument, with an optional second\n")
    sys.stderr.write("argument of the mode (add, drop, or full).\n\n")
    sys.stderr.write("Read the start of the script for more details.\n")
    sys.exit(1)

try:
    from Bio import SeqIO
except ImportError:
    sys.stderr.write("Required Biopython\n")
    sys.exit(1)


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


def add_or_remove_equals(ref_seq, read_seq, pos, cigar, add=True, drop=False):
    """Returns read_seq with equals signs for matched bases.

    Assumes both ref_seq and read_seq are using the same case.
    """
    if pos >= len(ref_seq):
        raise ValueError("Bad POS %i for a reference of length %i"
                         % (pos, len(ref_seq)))
    pending = read_seq
    answer = ""
    cigar_ops = decode_cigar(cigar)
    while cigar_ops:
        op_len, op = cigar_ops.pop(0)
        if op == "H":
            pass
        elif op == "S":
            answer += pending[:op_len]
            pending = pending[op_len:]
        elif op in "MX=":
            # print("%i%s - %s vs %s"
            #       % (op_len, op, ref_seq[pos:pos+op_len],pending[:op_len]))
            if op_len > len(pending):
                raise RuntimeError("Only %i bases left for %i%s" %
                                   (len(pending), op_len, op))
            if pos + op_len - 1 >= len(ref_seq):
                # TODO - Treat this as an error and terminate?
                sys.stderr.write("Warning, ran off end of %i bp reference by %i bp, pos %i, CIGAR %s\n"
                                 % (len(ref_seq), pos + op_len - len(ref_seq), pos, cigar))
                # e.g. NA12878.chromMT.SLX.maq.SRP000032.2009_07.bam read
                # SRR010937.2897481
                ref_seq += "?" * op_len  # Hack
            # for ref_base, read_base in
            # zip(ref_seq[pos:pos+op_len],pending[:op_len]):
            for i in range(op_len):
                read_base = pending[i]
                ref_base = ref_seq[pos + i]
                if ref_base == read_base or read_base == "=":
                    assert op != "X", "Bad CIGAR %s for %s" % (cigar, read_seq)
                    if add:
                        answer += "="
                    else:
                        # Remove any equals sign,
                        answer += ref_base
                else:
                    assert op != "=", "Bad CIGAR %s for %s" % (cigar, read_seq)
                    answer += read_base
            pending = pending[op_len:]
            pos += op_len
        elif op == "I":
            answer += pending[:op_len]
            pending = pending[op_len:]
        elif op in "DN":
            pos += op_len
        else:
            raise ValueError(
                "Unsupported CIGAR operator %s in %s" % (op, cigar))
    if pending:
        raise RuntimeError("Still had %i bases left" % len(pending))
    assert len(answer) == len(
        read_seq), "%s -> %s with %s" % (read_seq, answer, cigar)
    if drop and answer == "=" * len(answer) and len(answer) > 1:
        # Don't need the sequence at all!
        return "*", cigar.replace("M", "=")
    return answer, cigar

temp = add_or_remove_equals("ACGTWWWACGT", "TWW==", 3, "5M", True)
assert temp == ("=====", "5M"), temp
temp = add_or_remove_equals("ACGTWWWACGT", "TWW==", 3, "5M", True, True)
assert temp == ("*", "5="), temp
temp = add_or_remove_equals("ACGTWWWACGT", "TWW==", 3, "5=", True)
assert temp == ("=====", "5="), temp
temp = add_or_remove_equals("ACGTWWWACGT", "TWW==", 3, "5=", True, True)
assert temp == ("*", "5="), temp
temp = add_or_remove_equals("ACGTWWWACGT", "=====", 3, "5M", False)
assert temp == ("TWWWA", "5M"), temp
temp = add_or_remove_equals("ACGTWWWACGT", "=====", 3, "5=", False)
assert temp == ("TWWWA", "5="), temp
temp = add_or_remove_equals("ACGTWWWACGT", "==AWA", 3, "5M", True)
assert temp == ("==A==", "5M"), temp
temp = add_or_remove_equals("ACGTWWWACGT", "==AWA", 3, "2=1X2=", True)
assert temp == ("==A==", "2=1X2="), temp
temp = add_or_remove_equals("ACGTWWWACGT", "==A==", 3, "5M", False)
assert temp == ("TWAWA", "5M"), temp
temp = add_or_remove_equals("ACGTWWWACGT", "==A==", 3, "2=1X2=", False)
assert temp == ("TWAWA", "2=1X2="), temp

temp_mt = "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCTCCATGCATTTGGTATTTTCGTCTGGGGGGTGTGCACGCGATAGCATTGCGAGACGCTGGAGCCGGAGCACCCTATGTCGCAGTATCTGTCTTTGATTCCTGCCTCATTCTATTATTTATCGCACCTACGTTCAATATTACAGGCGAACATACCTACTAAAGTGTGTTAATTAATTAATGCTTGTAGGACATAATAATA"  # etc
temp = add_or_remove_equals(temp_mt, "GATCACAGGTCTATCACCCTATTAACCACTCACGGGAGCTCNCCAAT",
                            0, "47M", True)
assert ("=========================================N===AT", "47M") == temp, temp
temp = add_or_remove_equals(temp_mt, "TATTAACCCCTCACGTGATCTCTCCCTGCATTTTATTTTT",
                            19, "33M2D7M", True)
assert ("========C======T==T======C=============T", "33M2D7M") == temp, temp
# temp = add_or_remove_equals(temp_mt, "CGAAATCTGGTTCGTACTTCAGGGTCATAAAGCCTAAATAGCCCACCCGTTCCGCTTAGATAAGACATCACGATGG",
#                            86, "76M", True)
# assert "" == temp, temp
del temp_mt, temp


sys.stderr.write("Loading reference sequences from %s\n" % reference_filename)
try:
    import sqlite3
    del sqlite3
    reference = SeqIO.index_db(
        reference_filename + ".idx", reference_filename, "fasta")
except ImportError:
    reference = SeqIO.index(reference_filename, "fasta")
if not reference:
    sys.stderr.write(
        "No sequences found in FASTA reference file %s\n" % reference_filename)
    sys.exit(1)
sys.stderr.write("Sequences for %i reference(s) available\n" % len(reference))


ref_name = ""
ref_seq = ""
count = 0
mod = 0
bases = 0
for line in sys.stdin:
    if line[0] != "@":
        # Should be a read
        count += 1
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, rest = line.split(
            "\t", 11)
        if seq != "*":
            bases += len(seq)
            # TODO - Look at CIGAR or qual if SEQ is missing?
        if rname != "*" and not int(flag) & 0x4:
            # Mapped read
            if rname != ref_name:
                try:
                    ref_seq = str(reference[rname].seq).upper()
                except KeyError:
                    sys.stderr.write("Reference %s for read %s not in %s\n"
                                     % (rname, qname, reference_filename))
                    sys.exit(2)
                ref_name = rname
            # Add/remove equals signs in the read's sequence:
            try:
                seq, cigar = add_or_remove_equals(
                    ref_seq, seq, int(pos) - 1, cigar, add_equals, drop_seq)
                if seq == "*":
                    # According to spec, if omit SEQ must also omit QUAL (and
                    # samtools complains)
                    qual = "*"
            except Exception:
                sys.stderr.write(line)
                raise
            mod += 1
            line = "\t".join([qname, flag, rname, pos, mapq,
                              cigar, rnext, pnext, tlen, seq, qual, rest])
    sys.stdout.write(line)
sys.stderr.write("Modified %i out of %i reads\n" % (mod, count))
sys.stderr.write("In total %i bases in all %i reads\n" % (bases, count))
