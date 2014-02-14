#!/usr/bin/env python
"""Python script to convert SAM/BAM to SSPACE tab format.

Converts a SAM/BAM file containing mapped paired-end data into the
simple tab separated format used by the assembly scaffolder SSPACE:

    <contig1> <startpos_on_contig1> <endpos_on_contig1> <contig2> <startpos_on_contig2> <endpos_on_contig2>

e.g.

    contig1 100 150 contig1 350 300
    contig1 4000 4050 contig2 110 60

Essentially this is a rewrite of tools/sam_bam2Tab.pl script included
in SSPACE basic v2.0 to actually handle real SAM/BAM files where the
paired end data is correctly encoded using the FLAG field rather than
read name suffices.

Assuming your SAM/BAM file(s) have read groups, one tab file is created
for each read group - plus a libary file with the observed fragment
size information (taken from the TLEN field).

Simple usage with a paired-end SAM file:

$ ./sam_to_sspace_tab.py < original.sam converted

Simple usage with BAM files with conversion to/from SAM via samtools:

$ samtools view -h original.bam | ./sam_to_sspace_tab.py converted

This will produce files named converted_*.tab, one per read group
using the read group ID in the filename, plus converted.library
which is the main input file to give to SSPACE. Note the -h is
required with a BAM file in order to see the header information.

$ SSPACE_Basic_v2.0.pl -l converted.libraries -s original.fasta ...

TODO:

 * Autodetect orientation for library output
 * Actual mapped lengths (may need name sorted as in Perl original
   in order to efficiently get the partner read's mapped length)
 * Configurable size information in the output library file?
 * Report progress to stderr

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

def cigar_mapped_len(cigar):
    #TODO
    return 1

reads = 0
pairs = 0
interesting = 0
min_len = None
max_len = None
rg_handles = dict()
rg_lengths = dict()
for line in sys.stdin:
    if line[0]=="@":
        #Header line
        if line[1:3]=="RG":
            tags = line.rstrip().split("\t")
            rg = None
            for t in tags:
                if t.startswith("ID:"):
                    rg = t[3:]
            if rg is None:
                sys_exit("Missing ID in this read group line: %r" % line)
            rg_handles[rg] = open("%s_%s.tab" % (prefix, rg), "w")
            rg_lengths[rg] = dict()
        continue
    #Should be a read
    reads += 1
    qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, qual, tags = line.rstrip().split("\t", 11)
    rg_tags = [t for t in tags.split("\t") if t[:2] == "RG"]
    if not rg_tags:
        rg = None
        #Ignore this read? What about single library SAM/BAM files?
        continue
    elif len(rg_tags) > 1:
        sys_exit("Multiple RG tags in this line: %r" % line)
    else:
        rg = rg_tags[0]
        if not rg.startswith("RG:Z:"):
            sys_exit("Malformed RG tag %r in this line: %r" % (rg, line))
        rg = rg[5:]

    flag = int(flag)
    if (not (flag & 0x1) # Single end read
        or flag & 0x4 # Unmapped
        or flag & 0x8 # Partner unmapped
        or flag & 0x80 or not (flag & 0x40) #Only using read one (and the RNEXT/PNEXT information for read two)
        or flag & 0x100 or flag & 0x800 #Ignore secondary or supplementary alignments
        or flag & 0x200 # failed QC
        or flag & 0x400 # PCR or optical duplicate
        ):
        #Ignore this read
        continue
    len1 = cigar_mapped_len(cigar)
    len2 = len1 # TODO - this is a quick approximation...
    if flag & 0x16:
        # Read one is on the reverse strand
        end1 = int(pos)
        start1 = end1 - len1 + 1
    else:
        # Read one is on the forward strand
        start1 = int(pos)
        end1 = start1 + len1 - 1
    if flag & 0x20:
        # Partner (read two) is on the reverse strand
        end2 = int(pnext)
        start2 = end2 - len2 + 1
    else:
        # Partner (read two) is on the forward strand
        start2 = int(pnext)
        end2 = start2 + len2 - 1
    if rnext == "=":
        rnext = rname
    if rname == rnext:
        tlen = abs(int(tlen))
        if tlen:
            rg_lengths[rg][tlen] = rg_lengths[rg].get(tlen, 0) + 1
    else:
        interesting += 1

    try:
        handle = rg_handles[rg]
    except KeyError:
        sys_exit("Unexpected read group identifier %r in this line: %r" % (rg, line))

    handle.write("%s\t%i\t%i\t%s\t%i\t%i\n" % (rname, start1, end1, rnext, start2, end2))
    pairs += 1

for handle in rg_handles.values():
    handle.close()

sys.stderr.write("Extracted %i pairs from %i reads\n" % (pairs, reads))
sys.stderr.write("Of these, %i pairs are mapped to different contigs\n" % interesting)

handle = open(prefix + ".library", "w")
for rg in sorted(rg_lengths):
    lengths = rg_lengths[rg]
    size = 0
    error = 0.0
    direction = "FR" #TODO
    if lengths:
        print("Read group %s length range when mapped to same contig %i to %i, count %i, mean %0.1f"
              % (rg, min(lengths), max(lengths), len(lengths),
                 float(sum(lengths)) / len(lengths)))
        #This attempts to maximize pairings used (very inclusive)
        #TODO - Configurable?
        size = 0.5 * (min(lengths) + max(lengths))
        error = (max(lengths) - size) / size
        if error >= 1.0:
            #Ah. Can't cover all over them since SSPACE limits error to < 1.0 times size.
            size = float(sum(lengths)) / len(lengths) #median?
            error = 0.99999
    handle.write("%s TAB %s_%s.tab %i %0.5f %s\n" % (rg, prefix, rg, size, error, direction))
handle.close()
print("Now run SSPACE with your FASTA file and %s.library" % prefix)
