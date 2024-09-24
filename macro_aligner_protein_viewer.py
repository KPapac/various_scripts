#!/usr/bin/env python3
# -*- coding: utf-8 -*-
import os, io
import string
import numpy as np
import argparse

from Bio.Seq import Seq
from Bio.Align import MultipleSeqAlignment
from Bio import AlignIO, SeqIO

import panel as pn
import panel.widgets as pnw

pn.extension()

from bokeh.plotting import figure
from bokeh.models import ColumnDataSource, Plot, Grid, Range1d
from bokeh.models.glyphs import Text, Rect
from bokeh.layouts import gridplot
from bokeh.io import export_png

parser = argparse.ArgumentParser(
    prog='macro_aligner_viewer.py',
    description='Returns a png of a protein alignment.',
    epilog='Konstantinos Papachristos',
)
parser.add_argument('alignment', help='Path to the FASTA alignment file.')
args = parser.parse_args()


def get_colors(seqs):
    """make colors for bases in sequence"""
    text = [i for s in list(seqs) for i in s]
    # For protein sequences
    clrs = {
        'A': 'red',
        'C': 'blue',
        'D': 'green',
        'E': 'orange',
        'F': 'purple',
        'G': 'yellow',
        'H': 'pink',
        'I': 'brown',
        'K': 'cyan',
        'L': 'magenta',
        'M': 'lime',
        'N': 'olive',
        'P': 'gold',
        'Q': 'teal',
        'R': 'navy',
        'S': 'coral',
        'T': 'silver',
        'V': 'indigo',
        'W': 'maroon',
        'Y': 'salmon',
        'X': 'gray',  # Assuming gray for 'X' character
        '-': 'white',  # Assuming white for '-' character
    }
    colors = [clrs[i] for i in text]
    return colors


def view_alignment(aln, fontsize="9pt", plot_width=800):
    """Bokeh sequence alignment view"""
    # make sequence and id lists from the aln object
    seqs = [rec.seq for rec in (aln)]
    ids = [rec.id for rec in aln]
    text = [i for s in list(seqs) for i in s]
    colors = get_colors(seqs)
    N = len(seqs[0])
    S = len(seqs)
    width = .4
    x = np.arange(1, N + 1)
    y = np.arange(0, S, 1)
    # creates a 2D grid of coords from the 1D arrays
    xx, yy = np.meshgrid(x, y)
    # flattens the arrays
    gx = xx.ravel()
    gy = yy.flatten()
    # use recty for rect coords with an offset
    recty = gy + .5
    h = 1 / S
    # now we can create the ColumnDataSource with all the arrays
    source = ColumnDataSource(dict(x=gx, y=gy, recty=recty, text=text, colors=colors))
    plot_height = len(seqs) * 15 + 50
    x_range = Range1d(0, N + 1, bounds='auto')
    if N > 100:
        viewlen = 100
    else:
        viewlen = N
    # view_range is for the close up view
    view_range = (0, viewlen)
    tools = ""
    # entire sequence view (no text, with zoom)
    p = figure(
        title=f'{args.alignment}',
        plot_width=plot_width,
        plot_height=1500,
        x_range=(0, N),
        y_range=(0, S),
        tools=tools,
        min_border=5,
        toolbar_location='below',
    )
    rects = Rect(
        x="x",
        y="recty",
        width=1,
        height=1,
        fill_color="colors",
        line_color='black',
        line_width=0.1,
        fill_alpha=0.8,
    )
    p.add_glyph(source, rects)
    p.yaxis.visible = True
    p.grid.visible = True
    p = gridplot([[p]], toolbar_location='below')
    return p


aln = AlignIO.read(args.alignment, 'fasta')
p = view_alignment(aln, plot_width=1200)
export_png(p, filename=f'{os.path.basename(args.alignment)}.png')
