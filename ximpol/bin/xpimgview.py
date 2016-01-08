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


__description__ = 'Quick FITS image viewer'


from astropy.io import fits
import aplpy

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_no_grids
from ximpol.utils.logging_ import startmsg


def xpimgview(file_path):
    """Quick FITS image viewer.
    """
    hdu_list = fits.open(file_path)
    hdu_list.info()
    data = hdu_list[0].data
    with context_no_grids():
        fig = aplpy.FITSFigure(hdu_list[0], figure = plt.figure(0))
    fig.add_grid()
    fig.show_colorscale(cmap = 'afmhot')
    plt.show()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the input FITS image file')
    args = parser.parse_args()
    startmsg()
    xpimgview(args.file_path)
