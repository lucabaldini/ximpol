#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
# the Free Software Foundation; either version 3 of the License, or
# (at your option) any later version.
#
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
#
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.


"""matplotlib configuration module.
"""

import matplotlib
import numpy
import time
import os
from matplotlib import pyplot


from ximpol.__version__ import TAG
from ximpol import XIMPOL_DOC_FIGURES
from ximpol.utils.logging_ import logger


DEFAULT_FIG_WIDTH = 8.
DEFAULT_FIG_HEIGHT = 6.
DEFAULT_FIG_SIZE = (DEFAULT_FIG_WIDTH, DEFAULT_FIG_HEIGHT)


def property(key):
    """Return a given matplotlib configuration property.
    """
    return matplotlib.rcParams[key]

def context_two_by_two(scale=1.9):
    """Setup the current figure for a 2x2 panel.
    """
    _size = (scale*DEFAULT_FIG_WIDTH, scale*DEFAULT_FIG_HEIGHT)
    _rc = {'figure.figsize': _size}
    return matplotlib.rc_context(rc = _rc)

def context_no_grids():
    """Setup the current figure with no grids.
    """
    _rc = {'axes.grid': False}
    return matplotlib.rc_context(rc = _rc)

def overlay_tag(x=0.95, y=0.95, color='black'):
    """Overlay the ximpol tag on the current figure.
    """
    text = 'Created by ximpol %s on %s' % (TAG, time.asctime())
    pyplot.text(x, y, text, color=color, size=10, horizontalalignment='right',
                transform=pyplot.gca().transAxes)

def save_current_figure(file_name, folder=XIMPOL_DOC_FIGURES, clear=True,
                        show=False):
    """Save the current matplotlib figure in `XIMPOL_DOC_FIGURES`.

    Arguments
    ---------
    file_name : string
        The name of the output file.

    clear : bool
        If `True`, the current image is cleared after the fact.
    """
    file_path = os.path.join(folder, file_name)
    logger.info('Saving current figure to %s...' % file_path)
    pyplot.savefig(file_path, transparent=True)
    if show:
        pyplot.show()
    if clear:
        pyplot.clf()

def setup():
    """Basic setup.
    """
    matplotlib.rc('figure', facecolor='white', figsize=DEFAULT_FIG_SIZE)
    matplotlib.rc('lines', linewidth=2.0)
    matplotlib.rc('patch', linewidth=0.5, facecolor='blue', edgecolor='eeeeee',
                  antialiased=True)
    matplotlib.rc('text', hinting_factor=8, usetex=True)
    matplotlib.rc('mathtext', fontset='cm')
    matplotlib.rc('axes',facecolor='white', edgecolor='bcbcbc', grid=True,
                  titlesize='x-large', labelsize='large')
    matplotlib.rc('legend', fancybox=True, frameon=False, numpoints=1)



setup()
