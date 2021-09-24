#!/usr/bin/env python
"""Python script to render photos using bases A, C, G, and T for pixels.

Takes as input a PNG photo (JPEG should work if the right dependencies
are installed), and a FASTA sequence file, and uses them to produce a
PDF output image using ReportLab.

The motivation and example images are described on this blog post:
http://blastedbio.blogspot.co.uk/2013/08/pixelated-potato-posters-in-python.html
"""

from __future__ import print_function

import os

from Bio import SeqIO

from PIL import Image

import numpy as np

from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing, String
from reportlab.lib import colors
from reportlab.lib.units import cm, mm
from reportlab.pdfgen import canvas


# These are hueristic, and currently slighlty squash the image vertically,
# so h_scale should be a little smaller or the v_scale a little bigger.
h_scale = 0.140 * cm  # per bp
v_scale = 0.125 * cm  # per bp


def run(im, seq, pdf_file, main_caption):
    data = np.array(im)  # .reshape(shape_rgb, order="F")
    shape = data.shape[0:2]

    pixels = np.product(shape)
    print("Have %i base pairs, shape %r, and %i pixels" % (len(seq), shape, pixels))

    assert pixels <= len(seq)
    assert 0 <= data.min() <= data.max() <= 255

    # Open PDF
    width, height = page_size = h_scale * shape[1], v_scale * shape[0]
    print("Creating %s, %i by %i mm" % (pdf_file, width / mm, height / mm))
    c = canvas.Canvas(pdf_file, page_size)
    c.setTitle(main_caption)
    d = Drawing(*page_size)
    base = 0
    r_max = float(data[:, :, 0].max())
    g_max = float(data[:, :, 1].max())
    b_max = float(data[:, :, 2].max())
    for row in range(shape[0]):
        for col in range(shape[1]):
            # Ignore any alpha channel
            r, g, b = data[row, col, 0:3]
            # This scaling according to the channel maxima seemed like
            # a good idea to maximize contrast for use as a background
            # image:
            color = colors.Color(r / r_max, g / g_max, b / b_max)
            s = String(
                (col + 0.5) * h_scale,
                (shape[0] - row) * v_scale,
                seq[base],
                fillColor=color,
                fontSize=4,
                textAnchor="middle",
            )
            d.add(s)
            base += 1
    renderPDF.draw(d, c, 0, 0)
    c.showPage()
    c.save()


# print("A0: Suggest width %i, height %i pixels" % (841 * mm / h_scale, 1189 * mm / v_scale))
# --> A0: Suggest width 600, height 951 pixels
# print("A1: Suggest width %i, height %i pixels" % (594 * mm / h_scale, 841 * mm / v_scale))
# --> Suggest width 424, height 672 pixels
# print("A2: Suggest width %i, height %i pixels" % (420 * mm / h_scale, 594 * mm / v_scale))
# --> A2: Suggest width 300, height 475 pixels
# print("A3: Suggest width %i, height %i pixels" % (297 * mm / h_scale, 420 * mm / v_scale))
# --> A3: Suggest width 212, height 336 pixels
# print("A4: Suggest width %i, height %i pixels" % (210 * mm / h_scale, 297 * mm / v_scale))
# --> A4: Suggest width 150, height 237 pixels


# This is a hard coded list of potato images used for the
# poster backgrounds described here,
# http://blastedbio.blogspot.co.uk/2013/08/pixelated-potato-posters-in-python.html
for name, seq_file in [
    ("Purple on black 002", "chr06.fasta"),
    # ("Potato field 001", "chr07.fasta"),
    ("Potato field 002", "chr07.fasta"),
    ("Potato roots 001", "chr08.fasta"),
    ("Tractor 001", "chr09.fasta"),
    ("Potato flower", "chr01.fasta"),
    # ("Potato branch", "chr02.fasta"),
    ("Potato branch center", "chr02.fasta"),
    ("new branch 001", "chr02.fasta"),
    # ("Potato tubers", "chr03.fasta"),
    # ("Potato tubers2", "chr03.fasta"),
    ("Potato tubers 003", "chr03.fasta"),
    ("Potato leaves", "chr04.fasta"),
    ("Potato blue flowers", "chr05.fasta"),
    ("Blue Flower Brown", "chr05.fasta"),
    ("Blue Flower dark", "chr05.fasta"),
]:
    stem = name.lower().replace(" ", "_")
    png_file = "%s.png" % stem
    png_fileA = "%s_424x672.png" % stem
    png_fileB = "%s_600x951.png" % stem
    if not os.path.isfile(png_fileB):
        png_fileB = png_file
    if not os.path.isfile(png_fileA):
        png_fileA = png_fileB

    pdf_file = stem + "_%s.pdf"

    seq = str(SeqIO.read(seq_file, "fasta").seq)
    # Reduce runs of N to a single N to avoid visual distraction
    while "NNNNN" in seq:
        seq = seq.replace("NNNNN", "N")
    while "NN" in seq:
        seq = seq.replace("NN", "N")

    print("Drawing %s using %s" % (name, seq_file))
    for name, shape, png_file in [
        ("A4", (150, 237), png_fileB),
        ("A3", (212, 336), png_fileA),
        ("A2", (300, 475), png_fileB),
        ("A1", (424, 672), png_fileA),
        ("A0", (600, 951), png_fileB),
    ]:
        if not os.path.isfile(png_file):
            print("Missing %s" % png_file)
            continue
        if os.path.isfile(pdf_file % name):
            print("Skipping as %s exists..." % (pdf_file % name))
            continue
        print(
            "Size %s, using %i by %i pixels from %s"
            % (name, shape[0], shape[1], png_file)
        )
        im = Image.open(png_file).resize(shape)
        run(im, seq, pdf_file % name, name)
