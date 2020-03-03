"""
mpl_logoplot - logo.py

Plot logoplots in a matplotlib axis

"""

import os
import pickle

import numpy as np
import matplotlib as mpl

from . import paths


def logoplot(ax, seqmat, alphabet, height=None, color=None, x_offset=-0.5, y_offset=0.0, fontname='Ubuntu Mono'):
    font_paths = paths.get_paths(fontname, alphabet)

    if color is None:
        color = {}
    elif color in COLORS:
        color = COLORS[color]

    if height is None:
        height = np.ones(len(seqmat))

    y_upper = []
    y_lower = []

    patches = []
    minsize = 1e-4

    for i in range(len(seqmat)):
        o_pos = 0
        o_neg = 0

        patches.append([])

        for j in np.argsort(np.abs(seqmat[i])):
            letter = alphabet[j]
            verts, codes = font_paths[letter]
            verts = verts.copy()

            verts[:, 0] = verts[:, 0] + i + x_offset

            h = seqmat[i, j] * height[i]
            if h < 0:
                h = -h
                o_neg += h
                verts[:, 1] = verts[:, 1] * h - o_neg + y_offset
            else:
                verts[:, 1] = verts[:, 1] * h + o_pos + y_offset
                o_pos += h

            if h > minsize:
                path = mpl.path.Path(verts, codes)
                c = color.get(letter, '#666666')
                patch = mpl.patches.PathPatch(path, facecolor=c, lw=0)
                ax.add_patch(patch)

                patches[-1].append(patch)

        y_upper.append(o_pos)
        y_lower.append(o_neg)

    ax.set_xlim(x_offset, len(seqmat) + x_offset)

    ylow, _ = ax.set_ylim(-np.max(y_lower) + y_offset, np.max(y_upper) + y_offset)

    return patches


COLORS = {
    'dna': {
        'A': '#99cc33',
        'T': '#3366cc',
        'G': '#ff9900',
        'C': '#990000',
    },
    'protein': {
        'A': '#99cc33',
        'C': '#99cc33',
        'D': '#990000',
        'E': '#990000',
        'F': '#99cc33',
        'G': '#ff9900',
        'H': '#3366cc',
        'I': '#99cc33',
        'K': '#3366cc',
        'L': '#99cc33',
        'M': '#99cc33',
        'N': '#ff9900',
        'P': '#99cc33',
        'Q': '#ff9900',
        'R': '#3366cc',
        'S': '#ff9900',
        'T': '#ff9900',
        'V': '#99cc33',
        'W': '#99cc33',
        'Y': '#ff9900',
    }
}
