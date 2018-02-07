#!/usr/bin/env python
"""Python script for assembly comparison."""

from __future__ import print_function

import os
import sys
import warnings
from optparse import OptionParser

from Bio import BiopythonExperimentalWarning
with warnings.catch_warnings():
    warnings.simplefilter("ignore", BiopythonExperimentalWarning)
    from Bio import SearchIO

from Bio import SeqIO
from Bio.Blast.Applications import NcbiblastnCommandline
from Bio.Graphics import GenomeDiagram
from Bio.Graphics.GenomeDiagram import CrossLink
from Bio.SeqFeature import FeatureLocation, SeqFeature
from Bio.SeqIO.FastaIO import SimpleFastaParser

from reportlab.lib import colors
from reportlab.lib.units import cm

SPACER = 20000

MIN_GAP_JAGGY = 1000  # Sigils
MIN_GAP = 20000


usage = """Basic usage: multi_comparison.py assembly1.fasta assembly2.fasta assembly3.fasta -o figure.pdf

If a GenBank file exists next to FASTA file but with the extension *.gbk,
that will be loaded to draw any annotated genes. e.g. reference genome.

There should be a (nucleotide) BLAST database next to each FASTA
file, created with some thing like this such that the BLAST database files
are named reference.fasta.n* and the database is referenced simply as
reference.fasta when calling blastn:

$ makeblastdb -dbtype nucl -in assembly.fasta

The assembly FASTA files are expected to already be preordered, and
will each be drawn on one track with the contigs end to end (with a
spacer).
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
    elif pipe_name.startswith("gnl|") and pipe_name.count("|") == 2:
        return pipe_name.split("|")[2]
    else:
        return pipe_name


def sys_exit(msg, error_level=1):
    """Print error message to stdout and quit with given error level."""
    sys.stderr.write("%s\n" % msg.rstrip())
    sys.exit(error_level)

parser = OptionParser(usage=usage)
parser.add_option("-l", "--min-hit-len", dest="min_hit", type="int",
                  help="Minimum BLAST hit length to consider",
                  default=5000)
parser.add_option("-o", "--output", dest="pdf_filename",
                  help="Write PDF diagram to FILE (default automatic)",
                  default=None,
                  metavar="FILE")
(options, args) = parser.parse_args()
min_hit = int(options.min_hit)

if len(args) < 2:
    sys_exit("Requires two or more arguments!\n\n" + usage)
assemblies_fasta = args[:]
diagram_pdf = options.pdf_filename
if not diagram_pdf:
    sys_exit("Requires output PDF file to be specified!\n\n" + usage)
output_directory = os.path.dirname(diagram_pdf)

for assembly_fasta in assemblies_fasta:
    if not os.path.isfile(assembly_fasta):
        sys_exit("Assembly FASTA file not found: %r" % assembly_fasta)


def do_blast(query_fasta, db_fasta, blast_file):
    assert os.path.isfile(query_fasta)
    assert os.path.isfile(db_fasta)
    assert os.path.isfile(
        db_fasta + ".nhr"), "Missing database for %s" % db_fasta
    assert os.path.isfile(
        db_fasta + ".nin"), "Missing database for %s" % db_fasta
    assert os.path.isfile(
        db_fasta + ".nsq"), "Missing database for %s" % db_fasta
    cmd = NcbiblastnCommandline(query=query_fasta, db=db_fasta,
                                out=blast_file, outfmt=6,
                                evalue=1e-5)
    print(cmd)
    stdout, stderr = cmd()
    return


def filter_blast(blast_result, query_length):
    hsps = [hsp for hsp in blast_result.hsps if (
        hsp.query_end - hsp.query_start) >= min_hit]
    hsps = sorted(hsps, key=lambda hsp: hsp.hit_start)
    return blast_result.id, hsps


def add_jaggies(contig_seq, offset, gd_contig_features):
    """Add JAGGY features for any run of NNNN or XXXX in sequence."""
    contig_seq = contig_seq.upper().replace("X", "N")
    i = 0
    j = 0
    NNN = "N" * MIN_GAP_JAGGY
    while i < len(contig_seq):
        i = contig_seq.find(NNN, i)
        if i == -1:
            return
        j = i
        while j < len(contig_seq) and contig_seq[j] == "N":
            j += 1
        print("Adding jaggy")
        gd_contig_features.add_feature(SeqFeature(FeatureLocation(offset + i, offset + j)),
                                       sigil="JAGGY",
                                       color=colors.slategrey, border=colors.black)
        i = j + 1

max_len = 0
gd_diagram = GenomeDiagram.Diagram("Comparison")
reference_fasta = None
ref_offsets = dict()
gd_ref_features = None
for i, assembly_fasta in enumerate(assemblies_fasta):
    if not os.path.isfile(assembly_fasta):
        sys_exit("Assembly FASTA file not found: %r" % assembly_fasta)
    assembly_genbank = os.path.splitext(assembly_fasta)[0] + ".gbk"

    contig_offsets = dict()
    contig_lengths = dict()
    track_len = -SPACER
    with open(assembly_fasta) as h:
        for title, seq in SimpleFastaParser(h):
            contig_id = hack_ncbi_fasta_name(title.split(None, 1)[0])
            track_len += SPACER
            contig_offsets[contig_id] = track_len
            track_len += len(seq)
            contig_lengths[contig_id] = len(seq)
    # TODO - Configurable offset to visually align tracks?
    max_len = max(max_len, track_len)

    gd_track = gd_diagram.new_track(1 + 2 * i,
                                    name=assembly_fasta,
                                    greytrack=False, height=0.5,
                                    start=0, end=track_len)
    gd_contig_features = gd_track.new_set()

    blast_data = dict()
    if reference_fasta:
        output_stem = "%s_vs_%s" % (os.path.splitext(os.path.basename(assembly_fasta))[0],
                                    os.path.splitext(os.path.basename(reference_fasta))[0])
        output_stem = os.path.join(output_directory, output_stem)
        blast_file = output_stem + ".blast.tsv"

        if not os.path.isfile(blast_file):
            do_blast(assembly_fasta, reference_fasta, blast_file)
        print("Loading %s" % blast_file)
        for blast_result in SearchIO.parse(blast_file, "blast-tab"):
            blast_result.id = hack_ncbi_fasta_name(blast_result.id)
            contig_id, hits = filter_blast(
                blast_result, contig_lengths[blast_result.id])
            blast_data[contig_id] = hits
            # print("Using %i of %i hits for %s" % (len(hits), len(blast_result.hsps), contig_id))
    else:
        assert i == 0

    offset = 0
    if os.path.isfile(assembly_genbank):
        print("Drawing %s" % assembly_genbank)
        contigs = SeqIO.parse(assembly_genbank, "genbank")
    else:
        print("Drawing %s" % assembly_fasta)
        contigs = SeqIO.parse(assembly_fasta, "fasta")
    for contig in contigs:
        contig_id = hack_ncbi_fasta_name(contig.id)
        contig_len = len(contig)

        # Add feature for whole contig,
        loc = FeatureLocation(offset, offset + contig_len, strand=0)
        gd_contig_features.add_feature(SeqFeature(loc), color=colors.grey, border=colors.black,
                                       label=True, name=contig_id)
        # Mark any NNNN regions,
        add_jaggies(str(contig.seq), offset, gd_contig_features)
        # Mark any genes (if using GenBank file),
        for feature in contig.features:
            if feature.type != "gene":
                continue
            feature.location += offset
            gd_contig_features.add_feature(feature, sigil="BOX",
                                           color="lightblue", label=True,
                                           label_position="start",
                                           label_size=6, label_angle=0)

        # print "%s (len %i) offset %i" % (contig_id, contig_len, offset)

        if contig_id not in blast_data:
            # print("No BLAST matches for contig %s" % contig_id)
            offset += SPACER + contig_len
            continue
        blast_hsps = blast_data[contig_id]
        if not blast_hsps:
            # print("No BLAST matches for contig %s" % contig_id)
            offset += SPACER + contig_len
            continue
        # print("%i BLAST matches for contig %s" % (len(blast_hsps), contig_id))

        # Add cross-links,
        for hsp in blast_hsps:
            if hsp.hit_strand == -1:
                flip = True
                color = colors.blue
            else:
                flip = False
                color = colors.firebrick
            border = colors.lightgrey
            # Fade the colour based on percentage identify, 100% identity = 50%
            # transparent
            color = colors.Color(color.red, color.green,
                                 color.blue, alpha=(hsp.ident_pct / 200.0))
            assert offset == contig_offsets[hack_ncbi_fasta_name(hsp.query_id)]
            loc = FeatureLocation(offset + hsp.query_start,
                                  offset + hsp.query_end, strand=0)
            query = gd_contig_features.add_feature(
                SeqFeature(loc), color=color, border=border)
            try:
                r_offset = ref_offsets[hack_ncbi_fasta_name(hsp.hit_id)]
            except KeyError:
                if hack_ncbi_fasta_name(hsp.hit_id) != hsp.hit_id:
                    sys_exit("Could not find offset key %r for hit %r in dict (query id %r)"
                             % (hack_ncbi_fasta_name(hsp.hit_id), hsp.hit_id, hsp.query_id))
                else:
                    sys_exit("Could not find offset for hit %r in dict (query id %r)"
                             % (hsp.hit_id, hsp.query_id))
            loc = FeatureLocation(r_offset + hsp.hit_start,
                                  r_offset + hsp.hit_end, strand=0)
            hit = gd_ref_features.add_feature(
                SeqFeature(loc), color=color, border=border)
            gd_diagram.cross_track_links.append(
                CrossLink(query, hit, color, border, flip))

        offset += SPACER + contig_len

    # Ready for next pairwise comparison,
    reference_fasta = assembly_fasta
    ref_offsets = contig_offsets
    gd_ref_features = gd_contig_features

# Set size based on max track length?
page = (2 * cm + 5 * cm * len(assemblies_fasta), 100 * cm * max_len / 5000000)
gd_diagram.draw(format="linear", fragments=1,
                pagesize=page, start=0, end=max_len)
gd_diagram.write(diagram_pdf, "PDF")
print("Saved %r" % diagram_pdf)
