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

for contig in SeqIO.parse(assembly_fasta, "fasta"):
    if len(gd_diagram.tracks) > 10:
        #For testing, shortcut
        break
    contig_len = len(contig)
    gd_track_for_contig = gd_diagram.new_track(len(gd_diagram.tracks)+1,
                                               name=contig.id,
                                               greytrack=True, height=0.5,
                                               start=0, end=contig_len)
    gd_contig_features = gd_track_for_contig.new_set()
    for hit in blast_results[contig.id]:
        assert hit.query_id == contig.id
        assert hit.id == record.id or ("|%s|" % record.id) in hit.id, "%r vs %r" % (hit.id, record.id)
        for hsp in hit.hsps:
            color = colors.firebrick
            border = colors.lightgrey
            loc = FeatureLocation(hsp.query_start, hsp.query_end, strand=0)
            q = gd_contig_features.add_feature(SeqFeature(loc), color=color, border=border)
            loc = FeatureLocation(hsp.hit_start, hsp.hit_end, strand=0)
            h = gd_record_features.add_feature(SeqFeature(loc), color=color, border=border)
            gd_diagram.cross_track_links.append(CrossLink(q, h, color, border))
                                               

gd_feature_set = gd_track_for_features.new_set()
for feature in record.features:
    if feature.type != "gene":
        continue
    gd_feature_set.add_feature(feature, sigil="BIGARROW",
                               color="lightblue", label=True,
                               #name=str(i+1),
                               label_position="start",
                               label_size=6, label_angle=0)

page = ((3 + len(gd_diagram.tracks)) * cm, max_len * cm / 100000.0)
gd_diagram.draw(format="linear", pagesize=page, fragments=1,
                start=0, end=max_len)
gd_diagram.write(diagram_pdf, "PDF")
print "Saved %r" % diagram_pdf
