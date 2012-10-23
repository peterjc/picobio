#!/usr/bin/env python
"""Python script using a Bloom filter to select raw sequencing reads.

Often you can have sequence read data for a large genome (or even a
metagenome or other mixed sample), yet for a particular analysis be
interested only in reads mapping to a quite small reference set (e.g.
organelles like the mitochondria or chloroplast, or a symbiont, or
even just some particular genes of interest). You may want to try
several different read mapping protocols (e.g. a parameter sweep),
and doing this with the full read set would be too slow.

The idea of this tool is to act as a pre-filter, removing reads which
won't map, to concentrate only on those which might map.

The intention is to take as input a raw unaligned sequence read file
as input (in FASTQ, FASTA, SFF, or even unaligned SAM/BAM format),
and produce a filtered version as output.

You would also provide a set of linear and/or circular reference
sequences, the filtered output would only be reads which might map
to these reference sequences (the Bloom filter is probabilitistic,
but also we're only going to search for k-mers within each read,
not perform a full alignment).

We search reads against both the forward and reverse strand of the
reference. Currently this is done when building the filter, which
is problably faster at the cost of more memory (bigger filter). The
alternative would be to make a filter from the reference forward
strands only, and to check the reads and the reverse complement of
each read against it. I expect that to be slower.

TODO:

* Technically SFF support is easy via Biopython, but simple
  k-mer matching will be hampered by homopolymer errors.
  Also don't have to worry about handling paired reads.
* Command line switch to turn on paired read filtering
  (where if either read matches, both are kept in the output)
* Flexible FASTQ paired support to cope with paired or single
  (would negate need for another command line switch)
* FASTA paired support

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

try:
    from Bio.Seq import reverse_complement
except ImportError:
    sys_exit("Missing 'Bio' module, available from http://biopython.org")

try:
    import pydablooms
except ImportError:
    sys_exit("Missing 'dablooms' Python bindings, available from "
             "https://github.com/bitly/dablooms")

VERSION = "0.0.4"

def fasta_iterator(handle):
    """FASTA parser yielding (upper case sequence, raw record) string tuples."""
    raw = []
    seq = []
    for line in handle:
        if line.startswith(">"):
            if seq:
                yield "".join(seq).upper(), "".join(raw)
            seq = []
            raw = [line]
        elif raw:
            seq.append(line.strip())
            raw.append(line)
        else:
            raise ValueError("Bad FASTA line %r" % line)
    if raw:
        yield "".join(seq).upper(), "".join(raw)
    raise StopIteration

def fastq_iterator(handle):
    """FASTQ parser yielding (upper case sequence, raw record) string tuples.

    Note: Any text on the '+' line is dropped, so it isn't the exact raw
    record, it could be a slightly smaller file.
    """
    #TODO - Test this with nasty FASTQ files
    #Profile against reusing Biopython's FastqGeneralIterator
    while True:
        title = handle.readline()
        if not title:
            raise StopIteration
        if not title[0] == "@":
            raise ValueError("Expected FASTQ @ line, got %r" % title)
        seq = handle.readline()
        plus = handle.readline()
        if not plus[0] == "+":
            raise ValueError("Expected FASTQ + line, got %r" % plus)
        qual = handle.readline()
        if len(seq) != len(qual): #both include newline
            raise ValueError("Different FASTQ seq/qual lengths for %r" % title)
        yield seq.strip().upper(), title+seq+"+\n"+qual

def fastq_batched_iterator(handle):
    """FASTQ parser yielding (upper case sequence list, raw record(s) string) tuples.

    For use on interlaced paired FASTQ reads following the /1 and /2 suffix convention.
    """
    while True:
        #Read /1
        title = handle.readline()
        if not title:
            raise StopIteration
        if not title[0] == "@":
            raise ValueError("Expected FASTQ @ line, got %r" % title)
        id = title.split(None,1)[0]
        if not id.endswith("/1"):
            raise ValueError("Expected FASTQ record ending /1, got %r" % title)
        seq = handle.readline()
        plus = handle.readline()
        if not plus[0] == "+":
            raise ValueError("Expected FASTQ + line, got %r" % plus)
        qual = handle.readline()
        if len(seq) != len(qual): #both include newline
            raise ValueError("Different FASTQ seq/qual lengths for %r" % title)
        #Read /2
        title2 = handle.readline()
        if not title2[0] == "@":
            raise ValueError("Expected FASTQ @ line, got %r" % title2)
        id2 = title2.split(None,1)[0]
        if not id2.endswith("/2"):
            raise ValueError("Expected FASTQ record ending /2, got %r" % title2)
        if id[:-2] != id2[:-2]:
            raise ValueError("Expected paired FASTQ records, got %r and %r" % (title, title2))
        seq2 = handle.readline()
        plus = handle.readline()
        if not plus[0] == "+":
            raise ValueError("Expected FASTQ + line, got %r" % plus)
        qual2 = handle.readline()
        if len(seq) != len(qual): #both include newline
            raise ValueError("Different FASTQ seq/qual lengths for %r" % title2)        
        yield [seq.strip().upper(), seq2.strip().upper()], title+seq+"+\n"+qual+title2+seq2+"+\n"+qual2


def sam_iterator(handle):
    """SAM parser yielding (upper case sequence, raw record) string tuples.

    Checks reads are unmapped. Any header is discarded.
    """
    good_flags = set(["0", "77", "141"])
    for line in handle:
        if line[0]=="@":
            continue
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, rest = line.split("\t", 10)
        if flag in good_flags:
            yield seq.upper(), line
        else:
            sys_exit("Unexpected FLAG '%r' in SAM file, should be 0 (unmapped single read),\n"
                     "77 (0x4d, first of unmapped pair) or 141 (0x8d, second of unmapped pair).")

def sam_batched_iterator(handle):
    """SAM parser yielding (upper case sequence list, raw record(s) string) tuples.

    Checks reads are unmapped. Any header is discarded. Requires paired reads are
    consecutive in the file, FLAG 77 (0x4d) then FLAG 141 (0x8d).
    """
    good_flags = set(["0", "77", "141"])

    line = handle.readline()
    while line[0]=="@":
        line = handle.readline()

    while line:
        qname, flag, rname, pos, mapq, cigar, rnext, pnext, tlen, seq, rest = line.split("\t", 10)
        if flag == "0":
            #Unpaired unmapped read
            yield [seq.upper()], line
        elif flag == "77":
            #Paired read one
            line2 = handle.readline()
            qname2, flag2, rname2, pos2, mapq2, cigar2, rnext2, pnext2, tlen2, seq2, rest2 = line2.split("\t", 10)
            if qname != qname2:
                sys_exit("Missing second half of %s" % qname)
            if flag2 != "141":
                sys_exit("Expected FLAG 141 (0x8d) for second part of %s, got %r" % (qname, flag2))
            yield [seq.upper(), seq2.upper()], line+line2
        elif flag == "141":
            sys_exit("Missing first half of %s" % qname)
        else:
            sys_exit("Unexpected FLAG '%r' in SAM file, should be 0 (unmapped single read),\n"
                     "77 (0x4d, first of unmapped pair) or 141 (0x8d, second of unmapped pair).")
        line = handle.readline()
    raise StopIteration

def make_variants(seq, changes):
    #TODO - pass len(seq) == kmer as argument for speed?
    if changes > 1:
        raise NotImplementedError
    #Simple SNPs:
    for i in range(len(seq)):
        for letter in "ACGT":
            yield seq[:i] + letter + seq[i+1:]

def make_inserts(seq):
    """Given a k-mer, returns possible k-mers with single inserts."""
    for i in range(1,len(seq)):
        for letter in "ACGT":
            yield seq[:i] + letter + seq[i+1:-1]

def make_deletions(seq):
    """Given a (k+1)-mer, returns possible k-mers with single deletions."""
    for i in range(len(seq)):
        yield seq[:i] + seq[i+1:]

ambiguous_dna_values = {
    "A": "A",
    "C": "C",
    "G": "G",
    "T": "T",
    "M": "AC",
    "R": "AG",
    "W": "AT",
    "S": "CG",
    "Y": "CT",
    "K": "GT",
    "V": "ACG",
    "H": "ACT",
    "D": "AGT",
    "B": "CGT",
    "X": "GATC",
    "N": "GATC",
    }

def disambiguate(seq):
    good = set("ACGT")
    if good.issuperset(seq):
        yield seq
    else:
        assert len(seq) - sum(seq.count(letter) for letter in good) < 10, \
            "Have %i ambiguous chars in this %i sequence: %s" \
                % (len(seq) - sum(seq.count(letter) for letter in good), len(seq), seq)
        for i, letter in enumerate(seq):
            if letter not in good:
                for alt in ambiguous_dna_values[letter]:
                    #Now do some recursion...
                    #Would like to use Python 3.3+ syntax 'yield from ...' here
                    for new in disambiguate(seq[:i] + alt + seq[i+1:]):
                        yield new
                break

def build_filter(bloom_filename, linear_refs, circular_refs, kmer,
                 mismatches, inserts, deletions, error_rate=0.01, rc=True):
    #Using 5e-06 is close to a set for my example, both in run time
    #(a fraction more) and the number of reads kept (9528 vs 8058
    #with sets).
    simple = set()
    del_hashes = set()
    count = 0
    t0 = time.time()
    if linear_refs:
        for fasta in linear_refs:
            sys.stderr.write("Hashing linear references in %s\n" % fasta)
            handle = open(fasta)
            for upper_seq, raw_read in fasta_iterator(handle):
                #assert set(upper_seq).issubset("ACGT"), "%s contains %s" \
                #    % (raw_read.split("\n",1)[0], set(upper_seq).difference("ACGT"))
                #Note we do the disambiguate call on the fragments rather than
                #the whole reference to avoid too many levels of recursion.
                for i in range(0, len(upper_seq) - kmer):
                    for fragment in disambiguate(upper_seq[i:i+kmer]):
                        assert set(fragment).issubset("ACGT"), fragment
                        simple.add(fragment)
                        #bloom.add(fragment, kmer)
                        count += 1 #TODO - Can do this in one go from len(upper_seq)
                if deletions:
                    for i in range(0, len(upper_seq) - kmer+1):
                        for fragment in make_deletions(upper_seq[i:i+kmer+1]):
                            del_hashes.add(fragment)                    
            handle.close()

    if circular_refs:
        for fasta in circular_refs:
            sys.stderr.write("Hashing circular references in %s\n" % fasta)
            handle = open(fasta)
            for upper_seq, raw_read in fasta_iterator(handle):
                #assert set(upper_seq).issubset("ACGT"), "%s contains %s" \
                #    % (raw_read.split("\n",1)[0], set(upper_seq).difference("ACGT"))
                #Want to consider wrapping round the origin, add k-mer length:
                upper_seq += upper_seq[:kmer]
                for i in range(0, len(upper_seq) - kmer):
                    for fragment in disambiguate(upper_seq[i:i+kmer]):
                        assert set(fragment).issubset("ACGT"), fragment
                        simple.add(fragment)
                        #bloom.add(fragment, kmer)
                        count += 1 #TODO - Can do this in one go from len(upper_seq)
                if deletions:
                    for i in range(0, len(upper_seq) - kmer+1):
                        for fragment in make_deletions(upper_seq[i:i+kmer+1]):
                            del_hashes.add(fragment)
            handle.close()
    if rc:
        #Would popping be slow? Should mean less memory at once
        temp = simple.copy()
        for fragment in temp:
            simple.add(reverse_complement(fragment))
        del temp
    if mismatches or inserts or deletions:
        sys.stderr.write("Have %i unique k-mers before consider fuzzy matches\n" \
                         % (len(simple)))
        if deletions:
            #Do this first to avoid 3 large sets in memory!
            new = del_hashes
            del del_hashes
            new.update(simple)
            sys.stderr.write("Adding deletions brings this to %i unique k-mers\n" \
                             % len(new))
        else:
            new = simple.copy()
        if mismatches:
            for fragment in simple:
                for var in make_variants(fragment, mismatches):
                    new.add(var)
            sys.stderr.write("Adding %i mis-matches per k-mer, have %i unique k-mers\n" \
                             % (mismatches, len(new)))
        if inserts:
            for fragment in simple:
                for var in make_inserts(fragment):
                    new.add(var)
            sys.stderr.write("Adding inserts brings this to %i unique k-mers\n" \
                             % len(new))
        simple = new
    capacity = len(simple)
    bloom = pydablooms.Dablooms(capacity, error_rate, bloom_filename)
    for fragment in simple:
        bloom.add(fragment)
    bloom.flush()
    sys.stderr.write("Set and bloom filter of %i-mers created (%i k-mers considered, %i unique)\n" % (kmer, count, len(simple)))
    sys.stderr.write("Using Bloom filter with capacity %i and error rate %r\n" % (capacity, error_rate))
    sys.stderr.write("Building filters took %0.1fs\n" % (time.time() - t0))
    return simple, bloom

def go(input, output, format, paired, linear_refs, circular_refs, kmer, mismatches, inserts, deletions):
    if paired:
        if format=="fasta":
            #read_iterator = fasta_batched_iterator
            raise NotImplementedError
        elif format=="fastq":
            read_iterator = fastq_batched_iterator
        elif format=="sam":
            read_iterator = sam_batched_iterator
        else:
            sys_exit("Paired read format %r not recognised" % format)
    else:
        if format=="fasta":
            read_iterator = fasta_iterator
        elif format=="fastq":
            read_iterator = fastq_iterator
        elif format=="sam":
            read_iterator = sam_iterator
        else:
            sys_exit("Read format %r not recognised" % format)

    #Create new bloom file,
    handle, bloom_filename = tempfile.mkstemp(prefix="bloom-", suffix=".bin")
    #sys.stderr.write("Using %s\n" %  bloom_filename)
    simple, bloom = build_filter(bloom_filename, linear_refs, circular_refs,
                                 kmer, mismatches, inserts, deletions)

    #Now loop over the input, write the output
    if output:
        out_handle = open(output, "w")
    else:
        out_handle = sys.stdout

    if input:
        in_handle = open(input)
    else:
        in_handle = sys.stdin

    if format=="sam":
        out_handle.write("@HD\t1.4\tSO:unknown\n")

    in_count = 0
    out_count = 0
    t0 = time.time()
    filter_time = 0
    if paired:
        #If find a possible match in either of a pair of reads
        #keep them both (likewise for any multi-framgent set).
        for upper_seqs, raw_reads in read_iterator(in_handle):
            in_count += len(upper_seqs)
            wanted = False
            filter_t0 = time.time()
            for upper_seq in upper_seqs:
                for i in range(0, len(upper_seq) - kmer):
                    fragment = upper_seq[i:i+kmer]
                    if fragment in bloom and fragment in simple:
                        wanted = True
                        #Don't need to check rest of this read
                        break
                if wanted:
                    #Don't need to check the other reads
                    break
            filter_time += time.time() - filter_t0
            if wanted:
                out_handle.write(raw_reads)
                out_count += len(upper_seqs)
            if in_count % 100000 == 0:
                sys.stderr.write("Processed %i reads, kept %i (%0.1f%%), taken %0.1fs (of which %0.1fs in filter)\n" \
                                 % (in_count, out_count, (100.0*out_count)/in_count, time.time()-t0, filter_time))
    else:
        for upper_seq, raw_read in read_iterator(in_handle):
            in_count += 1
            wanted = False
            filter_t0 = time.time()
            for i in range(0, len(upper_seq) - kmer):
                fragment = upper_seq[i:i+kmer]
                #Can modify code to allow this syntax, see:
                #https://github.com/bitly/dablooms/pull/50
                #if bloom.check(fragment) and fragment in simple:
                if fragment in bloom and fragment in simple:
                    wanted = True
                    #Don't need to check rest of read
                    break
            filter_time += time.time() - filter_t0
            if wanted:
                out_handle.write(raw_read)
                out_count += 1
            if in_count% 100000 == 0:
                sys.stderr.write("Processed %i reads, kept %i (%0.1f%%), taken %0.1fs (of which %0.1fs in filter)\n" \
                                 % (in_count, out_count, (100.0*out_count)/in_count, time.time()-t0, filter_time))
    if input:
        in_handle.close()
    if output:
        out_handle.close()
    total_time = time.time() - t0
    sys.stderr.write("Running filter took %0.1fs, overhead %0.1fs, total %0.1fs\n" \
                         % (filter_time, total_time - filter_time, total_time))

    #Remove the bloom file
    del bloom
    os.remove(bloom_filename)
    sys.stderr.write("Kept %i out of %i reads (%0.1f%%)\n" % (out_count, in_count, out_count*100.0/in_count))

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
    #Matching
    parser.add_option("-k", "--kmer", dest="kmer",
                      type="int", metavar="KMER", default=35,
                      help="k-mer size for filtering (def. 35)")
    parser.add_option("-m", "--mismatches", dest="mismatches",
                      type="int", metavar="MM", default=0,
                      help="Number of mismatches per kmer (def. 0, max 1)")
    
    #Reads
    parser.add_option("-f", "--format", dest="format",
                      type="string", metavar="FORMAT", default="fasta",
                      help="Input (and output) read file format, one of 'fasta',"
                           " 'fastq' or 'sam' (unmapped reads only please).")
    #TODO - Make paired mode or single mode the default?
    parser.add_option("-i", "--input", dest="input_reads",
                      type="string", metavar="FILE",
                      help="Input file of unmapped reads to be filtered (def. stdin)")
    parser.add_option("-o","--output", dest="output_reads",
                      type="string", metavar="FILE",
                      help="Output file to write filtered reads to (def. stdout)")
    
    (options, args) = parser.parse_args()

    if len(sys.argv) == 1:
        parser.print_help()
        sys.exit(1)

    if not (10 <= options.kmer <= 100):
        parser.error("Using a k-mer value of %i is not sensible" % options.kmer)

    if options.mismatches < 0:
        parser.error("Number of mismatches per k-mer (here %i) cannot be negative" % options.mismatches)
    if options.mismatches > 1:
        parser.error("Number of mismatches per k-mer (here %i) currently limited to one" \
                     % options.mismatches)
    #TODO - Make substitions/inserts/deletions separate command line options?
    if options.mismatches:
        inserts = True
        deletions = True
    else:
        inserts = False
        deletions = False

    if (not options.linear_references) and (not options.circular_references):
        parser.error("You must supply some linear and/or circular references")

    if args:
        parser.error("No arguments expected")

    paired = True
    go(options.input_reads, options.output_reads, options.format, paired,
       options.linear_references, options.circular_references,
       options.kmer, options.mismatches, inserts, deletions)

if __name__ == "__main__":
    main()
