from Bio import SeqIO
import numpy as np
from PIL import Image

from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge
from reportlab.lib.units import mm, cm, inch
from reportlab.lib import colors

h_scale = 0.140 * cm #per bp
v_scale = 0.125 * cm #per bp

def run(im, seq, pdf_file, main_caption):
    data = np.array(im) #.reshape(shape_rgb, order="F") 
    shape = data.shape[0:2]

    pixels = np.product(shape)
    print "Have %i base pairs, shape %r, and %i pixels" % (len(seq), shape, pixels)
    
    assert pixels <= len(seq)
    assert 0 <= data.min() <= data.max() <= 255

    #Open PDF
    width, height = page_size = h_scale * shape[1], v_scale * shape[0]
    print "Creating %s, %i by %i mm" % (pdf_file, width / mm, height / mm)
    c = canvas.Canvas(pdf_file, page_size)
    c.setTitle(main_caption)
    d = Drawing(*page_size)
    base = 0
    r_max = float(data[:,:,0].max())
    g_max = float(data[:,:,1].max())
    b_max = float(data[:,:,2].max())
    for row in range(shape[0]):
        for col in range(shape[1]):
            #Ignore any alpha channel
            r, g, b = data[row, col, 0:3]
            color = colors.Color(r / r_max, g / g_max, b / b_max)
            s = String((col + 0.5) * h_scale, (shape[0]-row) * v_scale,
                       seq[base], fillColor = color,
                       fontSize = 4, textAnchor = "middle")
            d.add(s)
            base += 1
    renderPDF.draw(d, c, 0, 0)
    c.showPage()
    c.save()

"""
We want a portrait poster...

A0 841 x 1189 mm -> 600 x 951 pixels
A1 594 x 841 mm --> 424 x 672 pixels
A2 420 x 594 mm --> 300 x 475 pixels
A3 297 x 420 mm --> 212 x 336 pixels
A4 210 x 297 mm
"""

#print "A0: Suggest width %i, height %i pixels" % (841 * mm / h_scale, 1189 * mm / v_scale)
#--> A0: Suggest width 600, height 951 pixels
#print "A1: Suggest width %i, height %i pixels" % (594 * mm / h_scale, 841 * mm / v_scale)
#--> Suggest width 424, height 672 pixels
#print "A2: Suggest width %i, height %i pixels" % (420 * mm / h_scale, 594 * mm / v_scale)
#--> A2: Suggest width 300, height 475 pixels
#print "A3: Suggest width %i, height %i pixels" % (297 * mm / h_scale, 420 * mm / v_scale)
#--> A3: Suggest width 212, height 336 pixels

png_file = "potato_flower_424x672.png"
pdf_file = "potato_flower_%s.pdf"
main_caption = "Potato flower"

seq = str(SeqIO.read("chr01.fasta", "fasta").seq)
#Reduce runs of N to a single N to avoid visual distraction
while "NNNNN" in seq:
    seq = seq.replace("NNNNN", "N")
while "NN" in seq:
    seq = seq.replace("NN", "N")

for name, shape in [("A1", (424, 672)),
                    ("A2", (300, 475)),
                    ("A3", (212, 336))]:
    print "Size %s, using %i by %i pixels" % (name, shape[0], shape[1])
    im = Image.open(png_file).resize(shape)
    run(im, seq, pdf_file % name, main_caption)
