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


from matplotlib import pyplot

import matplotlib
import numpy


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

def setup():
    """Basic setup.
    """
    matplotlib.rc('figure', facecolor='white', figsize=DEFAULT_FIG_SIZE)
    matplotlib.rc('lines', linewidth=2.0)
    matplotlib.rc('patch', linewidth=0.5, facecolor='blue', edgecolor='eeeeee',
                  antialiased=True)
    matplotlib.rc('text', hinting_factor=8)
    matplotlib.rc('mathtext', fontset='cm')
    matplotlib.rc('axes',facecolor='white', edgecolor='bcbcbc', grid=True,
                  titlesize='x-large', labelsize='large')
    matplotlib.rc('legend', fancybox=True)


setup()
