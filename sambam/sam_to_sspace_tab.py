#!/usr/bin/env python
"""Python script to convert SAM/BAM to SSPACE tab format.

This is a *PROTOTYPE ONLY* which has had limited testing!

Converts a SAM/BAM file containing mapped paired-end data into the
simple tab separated format used by the assembly scaffolder SSPACE:

    <contig1> <startpos_on_contig1> <endpos_on_contig1> <contig2> <startpos_on_contig2> <endpos_on_contig2>

e.g.

    contig1 100 150 contig1 350 300
    contig1 4000 4050 contig2 110 60

Essentially this is a rewrite of tools/sam_bam2Tab.pl script included
in SSPACE basic v2.0 to actually handle real SAM/BAM files where the
paired end data is correctly encoded using the FLAG field rather than
read name suffices. It also generates a library file to use with this.

Assuming your SAM/BAM file(s) have read groups, one tab file is created
for each read group - plus a libary file with the observed fragment
size information (taken from the TLEN field).

Simple usage with a paired-end SAM file:

$ ./sam_to_sspace_tab.py < original.sam converted

Simple usage with BAM files with conversion to SAM via samtools:

$ samtools view -h original.bam | ./sam_to_sspace_tab.py converted

Note the -h is required with a BAM file in order to see the header
information.

This will produce files named converted_*.tab, one per read group
using the read group ID in the filename, plus converted.library
which is the main input file to give to SSPACE. This will attempt
to generate sensible values for the paired end insert size and
orientation, but you should check this and then run SSPACE:

$ SSPACE_Basic_v2.0.pl -l converted.libraries -s original.fasta ...

TODO:

 * Accept library information (size, orienation) via command line?
 * Output to a subdirectory? Would need relative paths...

Copyright Peter Cock 2014. All rights reserved. See:
https://github.com/peterjc/picobio
"""

import sys


def sys_exit(msg, error=1):
    sys.stderr.write(msg.rstrip() + "\n")
    sys.exit(error)


if len(sys.argv) != 2:
    sys_exit("Requires one argument, prefix for output tab files.")
prefix = sys.argv[1]


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


def cigar_mapped_len(cigar):
    """The aligned length is given by the sum of the CIGAR M/=/X/D/N operations."""
    if not cigar or cigar == "*":
        return 1  # Dummy value
    length = 0
    for op_length, op_code in decode_cigar(cigar):
        if op_code in "M=XDN":
            length += op_length
    return length


reads = 0
pairs = 0
interesting = 0
rg_handles = dict()
rg_lengths = dict()
rg_dir = dict()

cached = dict()  # Key by read name
for line in sys.stdin:
    if line[0] == "@":
        # Header line
        if line[1:3] == "RG":
            tags = line.rstrip().split("\t")
            rg = None
            for t in tags:
                if t.startswith("ID:"):
                    rg = t[3:]
            if rg is None:
                sys_exit("Missing ID in this read group line: %r" % line)
            rg_handles[rg] = open("%s_%s.tab" % (prefix, rg), "w")
            rg_lengths[rg] = []
            rg_dir[rg] = {"FR": 0, "RF": 0, "FF": 0}
        continue
    # Should be a read
    if reads % 500000 == 0:
        sys.stderr.write("Processed %i reads, %i pairs so far...\n" % (reads, pairs))
    reads += 1
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
        tags,
    ) = line.rstrip().split("\t", 11)
    rg_tags = [t for t in tags.split("\t") if t[:2] == "RG"]
    if not rg_tags:
        rg = None
        # Ignore this read? What about single library SAM/BAM files?
        continue
    elif len(rg_tags) > 1:
        sys_exit("Multiple RG tags in this line: %r" % line)
    else:
        rg = rg_tags[0]
        if not rg.startswith("RG:Z:"):
            sys_exit("Malformed RG tag %r in this line: %r" % (rg, line))
        rg = rg[5:]

    flag = int(flag)
    if (
        not (flag & 0x1)
        or flag & 0x4  # Single end read
        or flag & 0x8  # Unmapped
        or  # Partner unmapped
        # Neither R1 nor R2 (i.e. more than 2 parts)
        (flag & 0x40 and flag & 0x80)
        or not (flag & 0x40 or flag & 0x80)
        or flag & 0x100  # Unknown fragment number
        or flag & 0x800
        or flag & 0x200  # Ignore secondary or supplementary alignments
        or flag & 0x400  # failed QC
    ):  # PCR or optical duplicate
        # Ignore this read
        continue

    if rnext == "=":
        rnext = rname
    if qname in cached:
        # This is the second half of the pair (by file order)
        other_flag, other_rname, other_pos, other_cigar = cached.pop(qname)
        if other_rname != rnext or other_pos != pnext:
            sys.stderr.write(
                "Mapping position mismatch %s:%s versus %s:%s for %s\n"
                % (other_rname, other_pos, rnext, pnext, qname)
            )
            sys.stderr.write(line)
            sys_exit("Try running samtools fixmates?")
        if bool(flag & 0x10) != bool(other_flag & 0x20) or bool(flag & 0x20) != bool(
            other_flag & 0x10
        ):
            sys.stderr.write(
                "FLAG strand mismatch %i versus %i for %s\n" % (other_flag, flag, qname)
            )
            sys.stderr.write(line)
            sys_exit("Try running samtools fixmates?")
    else:
        # This is the first half of the pair (by file order), cache it
        cached[qname] = flag, rname, pos, cigar
        continue

    if flag & 0x40:
        # This is R1, other is R2
        assert other_flag & 0x80
        pass
    elif flag & 0x80:
        # This is R2, other is R1
        assert other_flag & 0x40
    else:
        assert False, "Bad FLAGs for %s (%i and %i)" % (qname, flag, other_flag)

    len1 = cigar_mapped_len(cigar)
    len2 = cigar_mapped_len(other_cigar)
    if flag & 0x10:
        # Read is on the reverse strand
        end1 = int(pos)
        start1 = end1 + len1 + 1
    else:
        # Read is on the forward strand
        start1 = int(pos)
        end1 = start1 + len1 - 1
    if flag & 0x20:
        # Partner (other read) is on the reverse strand
        end2 = int(pnext)
        start2 = end2 + len2 + 1
    else:
        # Partner (other read) is on the forward strand
        start2 = int(pnext)
        end2 = start2 + len2 - 1
    if rname == rnext:
        tlen = abs(int(tlen))
        if tlen:
            rg_lengths[rg].append(tlen)
            if (flag & 0x10) and (flag & 0x20):
                rg_dir[rg]["FF"] += 1
            elif flag & 0x10:
                assert end1 <= start1
                assert start2 <= end2
                # These are 'innies' --> <--,
                # Self:          end1 <---- start1
                # Other: start2 ----> end2
                #
                # Also consider overlapping reads as 'innies' --> <--
                # Self:      end1 <---- start1
                # Other: start2 ----> end2
                #
                # But these are 'outies' <-- -->
                # Self:  end1 <---- start1
                # Other:    start2 ----> end2
                if start1 < start2:
                    rg_dir[rg]["RF"] += 1  # 'outies' <-- -->
                else:
                    rg_dir[rg]["FR"] += 1  # 'innies' --> <--
            elif flag & 0x20:
                assert start1 <= end1
                assert end2 <= start2
                # Likewise, these are 'outies' <-- -->
                # Self:      start1 ----> end1
                # Other:  end2 <---- start2
                if start2 < start1:
                    rg_dir[rg]["RF"] += 1  # 'outies' <-- -->
                else:
                    rg_dir[rg]["FR"] += 1  # 'innies' --> <--
            else:
                rg_dir[rg]["FF"] += 1
    else:
        interesting += 1

    try:
        handle = rg_handles[rg]
    except KeyError:
        sys_exit("Unexpected read group identifier %r in this line: %r" % (rg, line))

    handle.write(
        "%s\t%i\t%i\t%s\t%i\t%i\n" % (rname, start1, end1, rnext, start2, end2)
    )
    pairs += 1

for handle in rg_handles.values():
    handle.close()

sys.stderr.write("Extracted %i pairs from %i reads\n" % (pairs, reads))
sys.stderr.write("Of these, %i pairs are mapped to different contigs\n" % interesting)
assert not cached, cached

handle = open(prefix + ".library", "w")
for rg in sorted(rg_lengths):
    lengths = rg_lengths[rg]
    size = 0
    error = 0.0
    direction = "??"
    if lengths:
        print(
            "Read group %s length range when mapped to same contig %i to %i, count %i, mean %0.1f"
            % (
                rg,
                min(lengths),
                max(lengths),
                len(lengths),
                float(sum(lengths)) / len(lengths),
            )
        )
        print(rg_dir[rg])
        assert sum(rg_dir[rg].values()) == len(lengths)
        # Pick most common direction
        direction = [
            d for d in rg_dir[rg] if rg_dir[rg][d] == max(rg_dir[rg].values())
        ][0]
        print("Most common pairing direction %s" % direction)
        # This attempts to maximize pairings used (very inclusive)
        # TODO - Configurable?
        size = 0.5 * (min(lengths) + max(lengths))
        error = (max(lengths) - size) / size
        if error >= 1.0:
            # Ah. Can't cover all over them since SSPACE limits error to < 1.0
            # times size.
            size = float(sum(lengths)) / len(lengths)  # median?
            error = 0.999
    handle.write(
        "%s TAB %s_%s.tab %i %0.3f %s\n" % (rg, prefix, rg, size, error, direction)
    )
handle.close()
print("Now run SSPACE with your FASTA file and %s.library" % prefix)
