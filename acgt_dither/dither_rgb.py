from Bio import SeqIO
import numpy as np
from PIL import Image

from reportlab.pdfgen import canvas
from reportlab.graphics import renderPDF
from reportlab.graphics.shapes import Drawing, String, Line, Rect, Wedge
from reportlab.lib.units import cm, inch
from reportlab.lib import colors

png_file = "potato_flower.png"
pdf_file = "potato_flower.pdf"
main_caption = "Potato flower"

#Load sequence
seq = str(SeqIO.read("SpaA1.fasta", "fasta").seq) * 10
shape = (400, 300)
shape_rgb = (400, 300, 4)
h_scale = 0.140 * cm #per bp
v_scale = 0.125 * cm #per bp

#Original is 1274 x 937 pixels, try about 20%
pixels = np.product(shape)
im = Image.open(png_file).resize(shape)
#im.show()
#data = im.getdata()
#assert len(data) == pixels, len(data)
assert shape == im.getbbox()[2:]
#data = np.array(data).reshape(shape, order="F")

data = np.array(im) #.reshape(shape_rgb, order="F") 
shape_rgb = data.shape
shape = shape_rgb[0:2]

assert shape_rgb == data.shape, "%r vs %r" % (shape, data.shape)
pixels = np.product(shape)
print "Have %i base pairs, shape %r, and %i pixels" % (len(seq), shape, pixels)

assert pixels <= len(seq)
assert 0 <= data.min() <= data.max() <= 255

#Open PDF
width, height = page_size = h_scale * shape[1], v_scale * shape[0]
c = canvas.Canvas(pdf_file, page_size)
c.setTitle(main_caption)
d = Drawing(*page_size)
base = 0
r_max = float(data[:,:,0].max())
g_max = float(data[:,:,1].max())
b_max = float(data[:,:,2].max())
for row in range(shape[0]):
    for col in range(shape[1]):
        r, g, b, alpha = data[row, col, :]
        #color = colors.Color((r - r_min) / (r_max - r_min),
        #                     (g - g_min) / (g_max - g_min),
        #                     (b - b_min) / (b_max - b_min))
        color = colors.Color(r / r_max, g / g_max, b / b_max)
        #color = colors.CMYKColor(black = (255 - data[col, row]) / 255.0)
        #From top left?
        s = String((col + 0.5) * h_scale, (shape[0]-row) * v_scale,
                   seq[base], fillColor = color,
                   fontSize = 4, textAnchor = "middle")
        d.add(s)
        base += 1
renderPDF.draw(d, c, 0, 0)
c.showPage()
c.save()
