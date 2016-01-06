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

def test():
    """
    """
    matplotlib.rc('lines', linewidth=2, color='r')
    matplotlib.rc('axes', linewidth=1.5, grid=True, labelsize='large')
    matplotlib.rc('grid', color='gray', linestyle=':', linewidth=0.5, alpha=1.0)
    matplotlib.rc('figure', facecolor='white')
    matplotlib.rc('text', usetex=False)

def bmh():
    """
    """
    matplotlib.rc('figure', facecolor='white')
    matplotlib.rc('lines', linewidth=2.0)
    matplotlib.rc('patch', linewidth=0.5, facecolor='blue', edgecolor='eeeeee',
                  antialiased=True)
    matplotlib.rc('text', hinting_factor=8)
    matplotlib.rc('mathtext', fontset='cm')
    matplotlib.rc('axes',facecolor='eeeeee', edgecolor='bcbcbc', grid=True,
                  titlesize='x-large', labelsize='large')
    matplotlib.rc('legend', fancybox=True)


bmh()
