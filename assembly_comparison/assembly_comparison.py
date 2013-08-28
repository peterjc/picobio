#!/usr/bin/env python
"""Python script for assembly comparison.
"""

import os
import sys

from Bio import SeqIO
from Bio import SearchIO
from Bio.SeqFeature import SeqFeature, FeatureLocation
from Bio.Blast.Applications import NcbiblastnCommandline

from reportlab.lib import colors
from reportlab.lib.colors import red, grey, orange, green, brown, blue, lightblue, purple
from reportlab.lib.units import cm

from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink

MIN_HIT = 5000
MIN_GAP = 50000

usage = """Usage: do_comparison.py assembly.fasta reference.fasta

If a reference GenBank file exists next to the reference FASTA file but
with the extension *.gbk, that will be loaded to draw any annotated genes.

There should be a (nucleotide) BLAST database next to the reference FASTA
file, created with some thing like this such that the BLAST database files
are named reference.fasta.n* and the database is referenced simply as
reference.fasta when calling blastn:

$ makeblastdb -in reference.fasta -dbtype nucl

"""

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

#TODO - Use argparse if API becomes non-trivial.
if len(sys.argv) != 3:
    stop_err(usage)
assembly_fasta, reference_fasta = sys.argv[1:]

reference_genbank = os.path.splitext(reference_fasta)[0] + ".gbk"
blast_file = assembly_fasta + ".blast.tsv"
diagram_pdf = assembly_fasta + ".blast.pdf"

if not os.path.isfile(assembly_fasta):
    stop_err("Assemlby FASTA file not found: %r" % assembly_fasta)

if not os.path.isfile(reference_fasta):
    stop_err("Reference FASTA file not found: %r" % reference_fasta)

if not os.path.isfile(reference_genbank):
    stop_err("Reference GenBank file not found: %r" % reference_genbank)

def do_blast(query_fasta, db_fasta, blast_file):
    assert os.path.isfile(query_fasta)
    assert os.path.isfile(db_fasta)
    assert os.path.isfile(db_fasta + ".nhr")
    assert os.path.isfile(db_fasta + ".nin")
    assert os.path.isfile(db_fasta + ".nsq")
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

#TODO - Multiple references (e.g. with plasmids)
if os.path.isfile(reference_genbank):
    record = SeqIO.read(reference_genbank, "genbank")
else:
    record = SeqIO.read(reference_fasta, "fasta")
max_len = len(record)

gd_diagram = GenomeDiagram.Diagram("Comparison")
gd_track_for_features = gd_diagram.new_track(1,
                                             name=record.name,
                                             greytrack=False, height=0.5,
                                             start=0, end=len(record))
gd_feature_set = gd_track_for_features.new_set()
#Add a dark grey background
gd_feature_set.add_feature(SeqFeature(FeatureLocation(0, len(record))),
                           sigil="BOX", color="grey", label=False),
for feature in record.features:
    if feature.type != "gene":
        continue
    gd_feature_set.add_feature(feature, sigil="BOX",
                               color="lightblue", label=True,
                               label_position="start",
                               label_size=6, label_angle=0)
gd_record_features = gd_track_for_features.new_set()


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
    hsps = [hsp for hsp in blast_result.hsps if (hsp.query_end - hsp.query_start) > MIN_HIT]
    hsps = sorted(hsps, key = lambda hsp: hsp.hit_start)
    plus = 0
    minus = 0
    flipped = False
    for hsp in hsps:
        if hsp.hit_strand == -1:
            minus += hsp.hit_end - hsp.query_start
        else:
            plus += hsp.hit_end - hsp.query_start
    if minus > plus:
        #Reverse the contig
        flipped = True
        hsps = [reverse_complement_hsp(hsp, query_length) for hsp in hsps]
        hsps = sorted(hsps, key = lambda hsp: hsp.hit_start)
    return make_offset(hsps), blast_result.id, hsps, flipped

def make_offset(blast_hsps):
    if not blast_hsps:
        return 0
    if blast_hsps[0].hit_strand == -1 and blast_hsps[-1].hit_strand == -1:
        #Assume whole thing flipped...
        offset = blast_hsps[0].hit_start - blast_hsps[-1].query_start
    else:
        offset = blast_hsps[0].hit_start - blast_hsps[0].query_start
    #return min(max(0, offset), max_len - contig_len)
    return offset

#Sort the contigs by horizontal position on the diagram
#(yes, this does mean parsing the entire BLAST output)
#(and yes, also the FASTA file to get the query lengths)
blast_data = sorted(filter_blast(b, len(contigs[b.id])) \
                        for b in SearchIO.parse(blast_file, "blast-tab"))


contig_tracks = []
for offset, contig_id, blast_hsps, flipped in blast_data:
    if not blast_hsps:
        continue

    #TODO - Use BLAST query length instead of parsing FASTA file?
    contig_len = len(contigs[contig_id])
    assert contig_len <= max_len
    offset = min(max(0, offset), max_len - contig_len)

    #Which track can we put this on?
    gd_track_for_contig = None
    gd_contig_features = None
    #print "%s needs offset %i" % (contig_id, offset)
    for track, fs in contig_tracks:
        #TODO - Can we calculate max of features instead of _used hack?
        if fs._used + MIN_GAP < offset:
            #Good, will fit in this track
            gd_track_for_contig = track
            gd_contig_features = fs
            break
    if not gd_track_for_contig:
        #print "Have %i tracks, adding one more" % len(contig_tracks)
        gd_track_for_contig = gd_diagram.new_track(2, #insert next to reference
                                                   name=contig_id,
                                                   greytrack=False, height=0.5,
                                                   start=0, end=max_len)
        gd_contig_features = gd_track_for_contig.new_set()
        contig_tracks.append((gd_track_for_contig, gd_contig_features))


    #Add feature for whole contig,
    loc = FeatureLocation(offset, offset + contig_len, strand=0)
    gd_contig_features.add_feature(SeqFeature(loc), color=colors.grey, border=colors.black,
                                   label=True, name=contig_id)
    gd_contig_features._used = offset +contig_len
    #print "%s (len %i) offset %i" % (contig_id, contig_len, offset)

    #Add cross-links,
    for hsp in blast_hsps:
        #print "%s:%i-%i hits %s:%i-%i" % (hsp.query_id, hsp.query_start, hsp.query_end,
        #                                  hsp.hit_id, hsp.hit_start, hsp.hit_end)
        if flipped:
            if hsp.hit_strand == -1:
                flip = True
                color = colors.darkgreen
            else:
                flip = False
                color = colors.purple
        else:
            if hsp.hit_strand == -1:
                flip = True
                color = colors.blue
            else:
                flip = False
                color = colors.firebrick
        border = colors.lightgrey
        #Fade the colour based on percentage identify, 100% identity = 50% transparent
        color = colors.Color(color.red, color.green, color.blue, alpha=(hsp.ident_pct/200.0))
        loc = FeatureLocation(offset + hsp.query_start, offset + hsp.query_end, strand=0)
        q = gd_contig_features.add_feature(SeqFeature(loc), color=color, border=border)
        loc = FeatureLocation(hsp.hit_start, hsp.hit_end, strand=0)
        h = gd_record_features.add_feature(SeqFeature(loc), color=color, border=border)
        gd_diagram.cross_track_links.append(CrossLink(q, h, color, border, flip))


page = (100*cm, 100*cm)
gd_diagram.draw(format="circular", circular=True, circle_core=0.5,
                pagesize=page, start=0, end=max_len)
gd_diagram.write(diagram_pdf, "PDF")
print "Saved %r" % diagram_pdf
