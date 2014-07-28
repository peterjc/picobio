#!/usr/bin/env python
"""Python script for assembly comparison.
"""

import os
import sys
import warnings
from optparse import OptionParser

from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import SearchIO

from Bio import SeqIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline

#TODO - make into command line options
MIN_HIT = 5000

usage = """Basic usage: order_assembly.py assembly.fasta reference.fasta output.fasta

There should be a (nucleotide) BLAST database next to the reference FASTA
file, created with some thing like this such that the BLAST database files
are named reference.fasta.n* and the database is referenced simply as
reference.fasta when calling blastn:

$ makeblastdb -dbtype nucl -in reference.fasta

Produces output.fasta based on the contigs in assembly.fasta which will
be reordered and/or reverse complemented to best match reference.fasta

Small repeat contigs which may match the reference in multiple locations
will hopefully be placed near one of these positions (typically the first).
"""

def hack_ncbi_fasta_name(pipe_name):
    """Turn 'gi|445210138|gb|CP003959.1|' into 'CP003959.1' etc.

    For use with NCBI provided FASTA and GenBank files to ensure
    contig names match up.

    Or Prokka's *.fna and *.gbk files, turning 'gnl|Prokka|contig000001'
    into 'contig000001'
    """
    if pipe_name.startswith("gi|") and pipe_name.endswith("|"):
        return pipe_name.split("|")[3]
    elif pipe_name.startswith("gnl|") and pipe_name.count("|")==2:
        return pipe_name.split("|")[2]
    else:
        return pipe_name

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

parser = OptionParser(usage=usage)
parser.add_option("-b", "--blast", dest="blast_filename",
                  help="Use/write BLAST tabular output to FILE (default automatic)",
                  default=None,
                  metavar="FILE")
parser.add_option("-m", "--min-len", dest="min_len", type="int",
                  help="Minimum contig length for FASTA output (if no BLAST hit)",
                  default=0)
#parser.add_option("-u", "--unmapped", dest="unmapped",
#                  help="Included unmapped contigs at the end",
#                  action="store_true")
(options, args) = parser.parse_args()

if len(args) != 3:
    stop_err("Requires three arguments!\n\n" + usage)
assembly_fasta, reference_fasta, output_fasta = args
blast_file = options.blast_filename
min_len = int(options.min_len)

output_stem = "%s_vs_%s" % (os.path.splitext(os.path.basename(assembly_fasta))[0],
                            os.path.splitext(os.path.basename(reference_fasta))[0])
output_stem = os.path.join(os.path.dirname(output_fasta), output_stem)

if not blast_file:
    blast_file = output_stem + ".blast.tsv"

if not os.path.isfile(assembly_fasta):
    stop_err("Assembly FASTA file not found: %r" % assembly_fasta)

if not os.path.isfile(reference_fasta):
    stop_err("Reference FASTA file not found: %r" % reference_fasta)

def do_blast(query_fasta, db_fasta, blast_file):
    assert os.path.isfile(query_fasta)
    assert os.path.isfile(db_fasta)
    if not (os.path.isfile(db_fasta + ".nhr") and \
            os.path.isfile(db_fasta + ".nin") and \
            os.path.isfile(db_fasta + ".nsq")):
        stop_err("Missing BLAST database for %s" % db_fasta)
    cmd = NcbiblastnCommandline(query=query_fasta, db=db_fasta,
                                out=blast_file, outfmt=6,
                                evalue=1e-5)
    print cmd
    stdout, stderr = cmd()
    return

if not os.path.isfile(blast_file):
    do_blast(assembly_fasta, reference_fasta, blast_file)

contigs = SeqIO.index(assembly_fasta, "fasta")
blast_results = SearchIO.index(blast_file, "blast-tab")

reference_parser = SeqIO.parse(reference_fasta, "fasta")

fasta_handle = open(output_fasta, "w")
fasta_saved_count = 0
fasta_short_dropped = 0

offset = 0
ref_offsets = dict()
for record in reference_parser:
    ref_offsets[hack_ncbi_fasta_name(record.id)] = offset
    offset += len(record)

def reverse_complement_hsp_fragment(frag, query_length):
    rev = SearchIO.HSPFragment(hit_id=frag.hit_id, query_id=frag.query_id)
    rev.query_start = query_length - frag.query_end
    rev.query_end = query_length - frag.query_start
    rev.hit_start = frag.hit_start
    rev.hit_end = frag.hit_end
    if frag.hit_strand == -1:
        rev.hit_strand = +1
    elif frag.hit_strand == +1:
        rev.hit_strand = -1
    else:
        #O or None,
        rev.hit_strand = frag.hit_strand
    return rev

def reverse_complement_hsp(hsp, query_length):
    rev = SearchIO.HSP(fragments = [reverse_complement_hsp_fragment(frag, query_length) \
                                    for frag in hsp.fragments[::-1]])
    rev.ident_pct = hsp.ident_pct
    return rev

def filter_blast(blast_result, query_length):
    hsps = [hsp for hsp in blast_result.hsps if (hsp.query_end - hsp.query_start) >= MIN_HIT]
    hsps = sorted(hsps, key = lambda hsp: hsp.hit_start)
    plus = 0
    minus = 0
    flipped = False
    for hsp in hsps:
        if hsp.hit_strand == -1:
            minus += hsp.hit_end - hsp.hit_start
        else:
            plus += hsp.hit_end - hsp.hit_start
    if minus > plus:
        #Reverse the contig
        flipped = True
        hsps = [reverse_complement_hsp(hsp, query_length) for hsp in hsps]
        hsps = sorted(hsps, key = lambda hsp: hsp.hit_start)
    return make_offset(hsps, query_length), blast_result.id, hsps, flipped


def weighted_median(values_and_weights):
    """Median of values with integer weights."""
    x = []
    count = sum(w for v, w in values_and_weights)
    map(x.extend,([v]*w for v, w in values_and_weights))
    return (x[count/2]+x[(count-1)/2])/2.


def make_offset(blast_hsps, contig_len):
    if not blast_hsps:
        return 0
    #Weighted by the HSP length:
    offset = int(weighted_median([(ref_offsets[hack_ncbi_fasta_name(hsp.hit_id)] + hsp.hit_start,
                                  hsp.hit_end - hsp.hit_start)
                                  for hsp in blast_hsps]))
    return offset


#Yes, this does end up parsing the entire FASTA file :(
contig_total_bp = sum(len(contigs[contig_id]) for contig_id in contigs)

#Sort the contigs by horizontal position on the diagram
#(yes, this does mean parsing the entire BLAST output)
#(and yes, also the FASTA file to get the query lengths)
blast_data = sorted(filter_blast(b, len(contigs[b.id])) \
                        for b in SearchIO.parse(blast_file, "blast-tab"))
contigs_shown = set()
contigs_shown_bp = 0
contig_tracks = []
for offset, contig_id, blast_hsps, flipped in blast_data:
    #TODO - Use BLAST query length instead of parsing FASTA file?
    contig = contigs[contig_id]
    contig_len = len(contig)
    if not blast_hsps:
        #Works, but only if contig appears in BLAST output at all
        #contigs_not_shown_bp += contig_len
        continue

    contigs_shown.add(contig_id)
    contigs_shown_bp += contig_len
    if contig_len < min_len:
        print("Note %s had BLAST hit but was only length %i" % (contig_id, contig_len))
    if flipped:
        SeqIO.write(contigs[contig_id].reverse_complement(id=True, name=True,
                                                          description="reversed"),
                    fasta_handle, "fasta")
    else:
        #Fast provided don't need to take reverse complement
        fasta_handle.write(contigs.get_raw(contig_id))
    fasta_saved_count += 1


#Now add the unmatched contigs
position = 0
unplaced = 0
for contig in SeqIO.parse(assembly_fasta, "fasta"):
    contig_id = contig.id
    if contig_id in contigs_shown:
        continue
    unplaced += 1
    contig_len = len(contig)
    if min_len <= contig_len:
        fasta_handle.write(contigs.get_raw(contig_id))
        fasta_saved_count += 1
    else:
        fasta_short_dropped += 1


assert unplaced == len(contigs) - len(contigs_shown), \
    "Only processed %i unplaced contigs, expected %i" % (unplaced, len(contigs) - len(contigs_shown))

print("Placed: %i of the %i contigs/scaffolds, %i bp"
      % (len(contigs_shown), len(contigs), contigs_shown_bp))
print("Unplaced: %i contigs/scaffolds, %i bp"
      % (len(contigs) - len(contigs_shown), contig_total_bp - contigs_shown_bp))
print("i.e. Placed %0.f%% of the assembly"
      % (contigs_shown_bp * 100.0 / contig_total_bp))

print("Wrote %i records to %r" % (fasta_saved_count, output_fasta))
print("Dropped %i short records" % fasta_short_dropped)
fasta_handle.close()
if fasta_saved_count + fasta_short_dropped != len(contigs):
    stop_err("Should have written %i records!" % (len(contigs) - fasta_short_dropped))
