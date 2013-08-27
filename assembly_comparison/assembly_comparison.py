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

def stop_err(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

#TODO - Use argparse if API becomes non-trivial.
if len(sys.argv) != 3:
    stop_err("Usage: do_comparison.py assembly.fasta reference.fasta")
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
record = SeqIO.read(reference_genbank, "genbank")
max_len = len(record)

gd_diagram = GenomeDiagram.Diagram("Comparison")
gd_track_for_features = gd_diagram.new_track(1,
                                             name=record.name,
                                             greytrack=True, height=1.0,
                                             start=0, end=len(record))
gd_record_features = gd_track_for_features.new_set()


def make_offset(blast_hsps):
    if blast_hsps[0].hit_strand == -1 and blast_hsps[-1].hit_strand == -1:
        #Assume whole thing flipped...
        offset = blast_hsps[0].hit_start - blast_hsps[-1].query_start
    else:
        offset = blast_hsps[0].hit_start - blast_hsps[0].query_start
    #return min(max(0, offset), max_len - contig_len)
    return offset

#Sort the contigs by horizontal position on the diagram
#(yes, this does mean parsing the entire BLAST output)
contig_offset_ids = sorted((make_offset(b.hsps), b.id)
                           for b in SearchIO.parse(blast_file, "blast-tab"))


contig_tracks = []
for offset, contig_id in contig_offset_ids:
    #TODO - Use BLAST query length instead of parsing FASTA file?
    contig_len = len(contigs[contig_id])
    assert contig_len <= max_len
    offset = min(max(0, offset), max_len - contig_len)

    blast_hits = blast_results[contig_id]
    assert len(blast_hits) == 1
    hit = blast_hits[0]
    assert hit.query_id == contig_id
    assert hit.id == record.id or ("|%s|" % record.id) in hit.id, "%r vs %r" % (hit.id, record.id)
    blast_hsps = blast_hits.hsps
    blast_hsps = [hsp for hsp in blast_hsps if (hsp.query_end - hsp.query_start) > MIN_HIT]
    if not blast_hsps:
        continue
    blast_hsps.sort(key = make_offset)
    del blast_hits, hit

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
    gd_contig_features.add_feature(SeqFeature(loc), color=colors.grey, border=colors.black)
    gd_contig_features._used = offset +contig_len
    #print "%s (len %i) offset %i" % (contig_id, contig_len, offset)

    #Add cross-links,
    for hsp in blast_hsps:
        #print "%s:%i-%i hits %s:%i-%i" % (hsp.query_id, hsp.query_start, hsp.query_end,
        #                                  hsp.hit_id, hsp.hit_start, hsp.hit_end)
        if hsp.hit_strand == -1:
            flip = True
            color = colors.blue
        else:
            flip = False
            color = colors.firebrick
        border = colors.lightgrey
        loc = FeatureLocation(offset + hsp.query_start, offset + hsp.query_end, strand=0)
        q = gd_contig_features.add_feature(SeqFeature(loc), color=color, border=border)
        loc = FeatureLocation(hsp.hit_start, hsp.hit_end, strand=0)
        h = gd_record_features.add_feature(SeqFeature(loc), color=color, border=border)
        gd_diagram.cross_track_links.append(CrossLink(q, h, color, border, flip))


gd_feature_set = gd_track_for_features.new_set()
for feature in record.features:
    if feature.type != "gene":
        continue
    gd_feature_set.add_feature(feature, sigil="BIGARROW",
                               color="lightblue", label=True,
                               #name=str(i+1),
                               label_position="start",
                               label_size=6, label_angle=0)

page = ((3 + 0.5*len(gd_diagram.tracks)) * cm, max_len * cm / 100000.0)
gd_diagram.draw(format="linear", pagesize=page, fragments=1,
                start=0, end=max_len)
gd_diagram.write(diagram_pdf, "PDF")
print "Saved %r" % diagram_pdf
