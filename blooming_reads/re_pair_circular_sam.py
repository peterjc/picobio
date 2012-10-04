#!/usr/bin/env python
"""Python script to re-pair sequencing reads to a circular genomes.

The idea is that you take your circular genomes, and double them,
and use a traditional mapper on that (in single end mode!), e.g.
I'm using mrfast at the moment for this. Then remove the duplicates
(mappings within the repeated region). Then run this script!

Input:

SAM file: Read name sorted SAM file, where the read names (QNAME)
have traditional FASTQ style /1 and /2 suffices indicating a
matched pair. Expects one SAM line where the mapping spans the
origin (with POS in range 1 to reference length, caculated
alignment end point would spill over).

FASTQ file: Used to add missing unmapped reads, often left out
in the SAM/BAM output of mapping tools. I want this to help
with downstream analysis.

Output:

SAM file: Still read name sorted, but with read pairings shown
via the FLAG field (the /1 and /2 are removed from the QNAME).
Additionally unmapped partners are present (from the FASTQ file
if not in the origin SAM file).

TODO:

Rename reads, moving /1 and /2 suffix into FLAG.

Adding missing unmapped reads.

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
from optparse import OptionParser

from Bio import SeqIO
from Bio.SeqIO.QualityIO import _get_sanger_quality_str as qual_str

def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

VERSION = "0.0.0"

solo0 = solo1 = solo2 = solo12 = 0


def go(input, output, raw_reads, linear_refs, circular_refs):

    if raw_reads:
        assert os.path.isfile(raw_reads)
        idx = raw_reads + ".idx"
        if os.path.isfile(idx):
            sys.stderr.write("Loading %s\n" % idx)
            raw = SeqIO.index_db(idx)
        else:
            sys.stderr.write("Creating %s\n" % idx)
            raw= SeqIO.index_db(idx, raw_reads, "fastq")
    else:
        raw = dict()
    sys.stderr.write("Have %i raw reads (used for unmapped partners)\n" % len(raw))
    
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
                if length == 2 * ref_len_circles[rname]:
                    #Return the length to its correct value (should have happened already)
                    sys.stderr.write("Fixing @SQ line for %s, length %i --> %i\n" % (rname, length, ref_len_circles[rname]))
                    line = "@SQ\tSN:%s\tLN:%i\n" % (rname, ref_len_circles[rname])
                else:
                    assert length == ref_len_circles[rname]
            elif rname is None:
                sys_exit("Bad @SQ line:\n%s" % line)
            else:
                sys_exit("This reference was not given!:\n%s" % line)
        output_handle.write(line)
        line = input_handle.readline()

    global solo0, solo1, solo2,solo12
    solo0 = solo1 = solo2 = solo12 = 0

    cur_read_name = None
    reads = set()
    while line:
        #SAM read
        qname, flag, rname, pos, rest = line.split("\t", 4)
        if " " in qname:
            #Stupid mrfast!
            qname = qname.split(None,1)[0]
        if rname in ref_len_circles and pos != "0":
            length = ref_len_circles[rname]
            if length <= int(pos) - 1:
                sys_exit("Have POS %s yet length of %s is %i (circular)\n" % (pos, rname, length))
        elif rname in ref_len_linear and pos != "0":
            length = ref_len_linear[rname]
            if length <= int(pos) - 1:
                        sys_exit("Have POS %s yet length of %s is %i (linear)\n" % (pos, rname, length))
        if qname[-2:] == "/1":
            qname = qname[:-2]
            frag = 1
        elif qname[-2:] == "/2":
            qname = qname[:-2]
            frag = 2
        else:
            frag = 0 #Assume unpaired
        if qname == cur_read_name:
            #Cache this, as a tuple - ordered to allow sorting on position:
            #Using a set will eliminate duplicates after adjusting POS
            reads.add((qname, frag, rname, pos, flag, rest))
        else:
            flush_cache(output_handle, reads, raw, ref_len_linear, ref_len_circles)
            reads = set([(qname, frag, rname, pos, flag, rest)])
            cur_read_name = qname
        #Next line...
        line = input_handle.readline()

    if reads:
        flush_cache(output_handle, reads, raw, ref_len_linear, ref_len_circles)

    if isinstance(input, basestring):
        input_handle.close()
    if isinstance(output, basestring):
        output_handle.close()

    sys.stderr.write("%i singletons; %i where only /1 mapped, %i where only /2 mapped, %i where both mapped\n" % (solo0, solo1, solo2, solo12))

def fixup_pairs(reads1, reads2, ref_len_linear, ref_len_circles):
    """Modify the two lists in-situ.

    TODO - Currently considers each reference in isolation!
    """
    if len(reads1) == len(reads2) == 1 and reads1[0][2] != reads2[0][2]:
        #Simple case, each maps to one place only, but they are on diff ref
        r1, r2 = make_mapped_pair(reads1[0], reads2[0])
        return [r1], [r2]
    fixed1 = []
    fixed2 = []
    refs1 = set(rname for qname, flag, rname, pos, rest in reads1)
    refs2 = set(rname for qname, flag, rname, pos, rest in reads2)
    for ref in sorted(refs1.union(refs2)):
        if ref in ref_len_linear:
            circular = False
            ref_lengths = ref_len_linear
        else:
            assert ref in ref_len_circles
            circular = True
            ref_lengths = ref_len_circles
        r1 = [(qname, flag, rname, pos, rest) for qname, flag, rname, pos, rest in reads1 if rname==ref]
        r2 = [(qname, flag, rname, pos, rest) for qname, flag, rname, pos, rest in reads2 if rname==ref]
        if ref in refs1 and ref in refs2:
            f1, f2 = fixup_same_ref_pairs(r1, r2, ref_lengths[ref], circular)
            fixed1.extend(f1)
            fixed2.extend(f2)
        else:
            #Mark them as mate unmapped (for now)
            #What if one other mapping only (so easy choice)?
            fixed1.extend((qname, flag | 0x41, rname, pos, rest) for qname, flag, rname, pos, rest in r1)
            fixed2.extend((qname, flag | 0x81, rname, pos, rest) for qname, flag, rname, pos, rest in r2)
    return fixed1, fixed2


def fixup_same_ref_pairs(reads1, reads2, length, circular):
    if len(reads1) == len(reads2) == 1:
        #Easy, can make them point at each other.
        #Need to look at locations & strands to decide if good mapping or not.
        r1, r2 = make_mapped_pair(reads1[0], reads2[0], happy=True)
        reads1 = [r1]
        reads2 = [r2]
    else:
        #TODO
        pass
    return reads1, reads2

def mark_mate(read_to_edit, mate_read, template_len=0, happy=False):
    qname1, flag1, rname1, pos1, rest1 = read_to_edit
    qname2, flag2, rname2, pos2, rest2 = mate_read

    assert qname1 == qname2

    if happy:
        #Set the properly paired bit
        flag1 |= 0x02

    mapq, cigar, rnext, pnext, tlen, etc = rest1.split("\t", 5)
    if rname1==rname2:
        rest1 = "\t".join([mapq, cigar, "=", pos2, str(template_len), etc])
    else:
        rest1 = "\t".join([mapq, cigar, rname2, pos2, str(template_len), etc])

    return [qname1, flag1, rname1, pos1, rest1]

def make_mapped_pair(read1, read2, template_len=0, happy=False):
    """Fill in RNEXT and PNEXT using each other's RNAME and POS."""
    read1 = mark_mate(read1, read2, template_len, happy)
    read2 = mark_mate(read2, read1, template_len, happy)
    return read1, read2

def flush_cache(handle, set_of_read_tuples, raw_dict, ref_len_linear, ref_len_circles):
    global solo0, solo1, solo2, solo12
    reads = sorted(set_of_read_tuples)
    if not reads:
        return

    #0x41 = 0x1 + 0x40 = paired, first in pair
    reads1 = [(qname, int(flag) | 0x41, rname, pos, rest) \
              for (qname, frag, rname, pos, flag, rest) \
              in reads if frag==1]

    #0x81 = 0x1 + 0x80 = paired, second in pair
    reads2 = [(qname, int(flag) | 0x81, rname, pos, rest) \
              for (qname, frag, rname, pos, flag, rest) \
              in reads if frag==2]

    reads0  = [(qname, int(flag), rname, pos, rest) \
              for (qname, frag, rname, pos, flag, rest) \
              in reads if frag==0]

    assert len(reads1) + len(reads2) + len(reads0) == len(reads)

    if reads0:
        assert not reads1 and not reads2, reads
        #All singletons
        solo0 += 1
    else:
        assert not reads0, reads
        assert reads1 or reads2, reads
        #Pairs
        if not reads1:
            solo2 += 1
            #0x8 = partner unmapped
            reads2 = [(qname, flag | 0x8, rname, pos, rest) \
                      for (qname, flag, rname, pos, rest) in reads2]
            #Assume first read2 is best one
            qname, rname, pos, flag, rest = reads2[0]
            flag = 0x1 + 0x4 + 0x40 #Paired, this is unmapped, first in pair
            rec = raw_dict[qname + "/1"]
            rest = "255\t*\t%s\t%s\t0\t%s\t%s\n" % (rname, pos, rec.seq, qual_str(rec))
            reads1 = [(qname, flag, "*", "0", rest)]
        elif not reads2:
            solo1 += 1
            #0x8 = partner unmapped
            reads1 = [(qname, flag | 0x8, rname, pos, rest) \
                      for (qname, flag, rname, pos, rest) in reads1]
            #Assume first read1 is best one:
            qname, rname, pos, flag, rest = reads1[0]
            flag = 0x1 + 0x4 + 0x80 #Paired, this is unmapped, second in pair
            rec = raw_dict[qname + "/2"]
            rest = "255\t*\t%s\t%s\t0\t%s\t%s\n" % (rname, pos, rec.seq, qual_str(rec))
            reads2 = [(qname, flag, "*", "0", rest)]
        else:
            solo12 += 1
            reads1, reads2 = fixup_pairs(reads1, reads2, ref_len_linear, ref_len_circles)

    for qname, flag, rname, pos, rest in reads0+reads1+reads2:
        #Assume that 'rest' has the trailing \n
        handle.write("\t".join([qname, str(flag), rname, pos, rest]))


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
    parser.add_option("-i", "--input", dest="input_reads",
                      type="string", metavar="FILE",
                      help="Input file of SAM format mapped reads to be processed (def. stdin)")
    parser.add_option("-r", "--reads", dest="raw_reads",
                      type="string", metavar="FILE",
                      help="Input file of FASTQ format unmapped reads (for finding unmapped partners)")
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
    go(options.input_reads, options.output_reads, options.raw_reads,
       options.linear_references, options.circular_references)

if __name__ == "__main__":
    main()
