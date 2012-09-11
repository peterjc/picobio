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

VERSION = "0.0.1"

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
    print "Linear references: %r" % options.linear_references
    print "Circular references: %r" % options.circular_references
    print "Using k-mer size %i" % options.kmer

if __name__ == "__main__":
    main()
