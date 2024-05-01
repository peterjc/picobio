from __future__ import print_function

import numpy as np
from Bio import SeqIO
from PIL import Image
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing
from reportlab.graphics.shapes import String
from reportlab.lib import colors
from reportlab.lib.units import cm
from reportlab.pdfgen import canvas

png_file = "Swanson_et_al_2012_fig1a.png"
pdf_file = "Swanson_et_al_2012_fig1a.pdf"
main_caption = "Swanson et al (2012) Figure 1"

# Load sequence
seq = SeqIO.read("SpaA1.fasta", "fasta").seq
shape = (239, 176)
scale = 0.125 * cm  # per bp

# Original is 1274 x 937 pixels, try about 20%
pixels = np.product(shape)
im = Image.open(png_file).resize(shape)
# im.show()
data = im.getdata()
assert len(data) == pixels, len(data)
assert shape == im.getbbox()[2:]
data = np.array(data).reshape(shape, order="F")
assert shape == data.shape
pixels = np.product(shape)
print("Have %i base pairs, and %i pixels" % (len(seq), pixels))

assert pixels <= len(seq)
assert 0 <= data.min() <= data.max() <= 255

# Open PDF
width, height = page_size = [x * scale for x in shape]
c = canvas.Canvas(pdf_file, page_size)
c.setTitle(main_caption)
d = Drawing(*page_size)
base = 0
for row in range(shape[1]):
    for col in range(shape[0]):
        color = colors.CMYKColor(black=(255 - data[col, row]) / 255.0)
        # From top left?
        s = String(
            (col + 0.5) * scale,
            (shape[1] - row) * scale,
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
