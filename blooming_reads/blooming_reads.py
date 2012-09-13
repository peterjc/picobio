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

The intension is to take as input a raw unaligned sequence read file
as input (in FASTQ, FASTA, SFF, or even unaligned SAM/BAM format),
and produce a filtered version as output.

You would also provide a set of linear and/or circular reference
sequences, the filtered output would only be reads which might map
to these reference sequences (the Bloom filter is probabilitistic,
but also we're only going to search for k-mers within each read,
not perform a full alignment).
"""
import sys
import os
import tempfile
import time
from optparse import OptionParser

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg)
    sys.exit(error_level)

try:
    import pydablooms
except ImportError:
    sys_exit("Missing 'dablooms' Python bindings, available from "
             "https://github.com/bitly/dablooms")

VERSION = "0.0.2"

def fasta_iterator(filename):
    """FASTA parser yielding (upper case sequence, raw record) string tuples."""
    if isinstance(filename, basestring):
        handle = open(filename)
    else:
        handle = filename
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
    if isinstance(filename, basestring):
        handle.close()
    if raw:
        yield "".join(seq).upper(), "".join(raw)
    raise StopIteration

def build_filter(bloom_filename, linear_refs, circular_refs, kmer,
                 error_rate=0.0005):
    #Using 5e-06 is close to a set for my example, both in run time
    #(a fraction more) and the number of reads kept (9528 vs 8058
    #with sets).
    simple = set()
    count = 0
    t0 = time.time()
    if linear_refs:
        for fasta in linear_refs:
            sys.stderr.write("Hashing linear references in %s\n" % fasta)
            for upper_seq, raw_read in fasta_iterator(fasta):
                for i in range(0, len(upper_seq) - kmer):
                    fragment = upper_seq[i:i+kmer]
                    simple.add(fragment)
                    #bloom.add(fragment, kmer)
                    count += 1 #TODO - Can do this in one go from len(upper_seq)

    if circular_refs:
        for fasta in circular_refs:
            sys.stderr.write("Hashing circular references in %s\n" % fasta)
            for upper_seq, raw_read in fasta_iterator(fasta):
                #Want to consider wrapping round the origin, add k-mer length:
                upper_seq += upper_seq[:kmer]
                for i in range(0, len(upper_seq) - kmer):
                    fragment = upper_seq[i:i+kmer]
                    simple.add(fragment)
                    #bloom.add(fragment, kmer)
                    count += 1 #TODO - Can do this in one go from len(upper_seq)

    capacity = len(simple)
    bloom = pydablooms.Dablooms(capacity, error_rate, bloom_filename)
    for fragment in simple:
        bloom.add(fragment)
    bloom.flush()
    sys.stderr.write("Set and bloom filter of %i-mers created (%i k-mers considered, %i unique)\n" % (kmer, count, len(simple)))
    sys.stderr.write("Using Bloom filter with capacity %i and error rate %r\n" % (capacity, error_rate))
    sys.stderr.write("Building filters took %0.1fs\n" % (time.time() - t0))
    return simple, bloom

def go(input, output, format, linear_refs, circular_refs, kmer):
    if format=="fasta":
        read_iterator = fasta_iterator
    elif format=="fastq":
        raise NotImplementedError
    elif format=="sam":
        raise NotImplementedError
    else:
        sys_exit("Read format %r not recognised" % format)

    #Create new bloom file,
    handle, bloom_filename = tempfile.mkstemp(prefix="bloom-", suffix=".bin")
    print bloom_filename
    simple, bloom = build_filter(bloom_filename, linear_refs, circular_refs, kmer)

    #Now loop over the input, write the output
    if output:
        out_handle = open(output, "w")
    else:
        out_handle = sys.stdout
    if not input:
        input = sys.stdin

    in_count = 0
    out_count = 0
    t0 = time.time()
    filter_time = 0
    for upper_seq, raw_read in read_iterator(input):
        in_count += 1
        wanted = False
        filter_t0 = time.time()
        for i in range(0, len(upper_seq) - kmer):
            fragment = upper_seq[i:i+kmer]
            #Can modify code to allow this syntax, see:
            #https://github.com/bitly/dablooms/pull/50
            #if bloom.check(fragment) and fragment in simple:
            if fragment in bloom: # and fragment in simple:
                wanted = True
                #Don't need to check rest of read
                break
        filter_time += time.time() - filter_t0
        if wanted:
            out_handle.write(raw_read)
            out_count += 1
    if output:
        out_handle.close()
    total_time = time.time() - t0
    sys.stderr.write("Running filter took %0.1fs, overhead %0.1fs, total %0.1fs\n" \
                         % (filter_time, total_time - filter_time, total_time))

    #Remove the bloom file
    del bloom
    os.remove(bloom_filename)
    sys.stderr.write("Kept %i out of %i reads\n" % (out_count, in_count))

def main():
    parser = OptionParser(usage="usage: %prog [options]",
                          version="%prog "+VERSION)
    parser.add_option("-l", "--lref", dest="linear_references",
                      type="string", metavar="FILE", action="append",
                      help="""FASTA file of linear reference sequence(s)
                           Several files can be given if required.""")
    parser.add_option("-c", "--cref", dest="circular_references",
                      type="string", metavar="FILE", action="append",
                      help="""FASTA file of circular reference sequence(s)
                           Several files can be given if required.""")
    parser.add_option("-k", "--kmer", dest="kmer",
                      type="int", metavar="KMER", default=35,
                      help="k-mer size for filtering (def. 35)")
    parser.add_option("-f", "--format", dest="format",
                      type="string", metavar="FORMAT", default="fasta",
                      help="Input (and output) read file format, one of 'fasta',"
                           " 'fastq' or 'sam' (unmapped reads only please).")
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

    if (not options.linear_references) and (not options.circular_references):
        parser.error("You must supply some linear and/or circular references")

    if args:
        parser.error("No arguments expected")

    print "Will attempt to filter %s --> %s" % (options.input_reads, options.output_reads)
    #print "Linear references: %r" % options.linear_references
    #print "Circular references: %r" % options.circular_references
    #print "Using k-mer size %i" % options.kmer
    go(options.input_reads, options.output_reads, options.format,
       options.linear_references, options.circular_references, options.kmer)

if __name__ == "__main__":
    main()
