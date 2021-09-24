#!/usr/bin/env python
usage = """Python script to restore SEQ entries recorded as just * in SAM/BAM.

This script is designed to be used on the BWA-MEM output where it seems
with the -a option additional alignments are reported with the SEQ and
QUAL fields set to just *.

If the alignment was recorded with SEQ * and a CIGAR string with hard
clipping, we restore the full SEQ and alter the CIGAR to say soft trimming
was used.

Developed with the output from:

    ~/Downloads/bwa-0.7.10/bwa mem -p -S -a REF.fas PAIRED.fastq > PAIRED.sam

and:

    ~/Downloads/bwa-0.7.10/bwa mem -S -a REF.fas SINGLE.fastq > SINGLE.sam

This script is designed to be used as part of a Unix pipeline. It reads
SAM format data from stdin, and writes SAM format data to stdout.

TODO: May need to refer to the original unaligned FASTQ file in some
cases (e.g. if the BWA output was sorted), but the assumption is that
the first mapped read will include the full SEQ and QUAL fields.

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys


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
            raise ValueError("Invalid character %s in CIGAR %s" % (letter, cigar))
    return answer


assert decode_cigar("14S15M1P1D3P54M1D34M5S") == [
    (14, "S"),
    (15, "M"),
    (1, "P"),
    (1, "D"),
    (3, "P"),
    (54, "M"),
    (1, "D"),
    (34, "M"),
    (5, "S"),
]


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


last_seq = None
last_qual = None
last_name = None
last_frag = 0
count = 0
mod = 0
for line in sys.stdin:
    if line[0] != "@":
        # Should be a read
        count += 1
        (
            qname,
            flag,
            rname,
            pos,
            mapq,
            cigar,
            rnext,
            pnext,
            tlen,
            seq,
            qual,
            rest,
        ) = line.split("\t", 11)
        if seq == "*":
            if qname == last_name and get_frag(flag) == last_frag:
                assert last_seq
                exp_len = cigar_seq_len(cigar)
                if exp_len < len(last_seq) and "H" in cigar:
                    # Ought to work if record it as sort trimming...
                    cigar = cigar.replace("H", "S")
                    exp_len = cigar_seq_len(cigar)
                assert exp_len == len(last_seq), (
                    "Cached SEQ %r length %i, but this read CIGAR expects length %i:\n%s"
                    % (last_seq, len(last_seq), cigar_seq_len(cigar), line)
                )
                seq = last_seq
                if qual == "*":
                    qual = last_qual
                mod += 1
                line = "\t".join(
                    [
                        qname,
                        flag,
                        rname,
                        pos,
                        mapq,
                        cigar,
                        rnext,
                        pnext,
                        tlen,
                        seq,
                        qual,
                        rest,
                    ]
                )
        elif "H" not in cigar:
            # Cache the SEQ
            last_name = qname
            last_frag = get_frag(flag)
            last_seq = seq
            last_qual = qual
    sys.stdout.write(line)
sys.stderr.write("Modified %i out of %i reads\n" % (mod, count))
