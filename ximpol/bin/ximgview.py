#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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



import matplotlib.pyplot as plt
from astropy.io import fits
import aplpy

from ximpol.__logging__ import logger, startmsg



def ximgview(filePath):
    """ Quick FITS image viewer.
    """
    hdulist = fits.open(filePath)
    hdulist.info()
    data = hdulist[0].data
    fig = aplpy.FITSFigure(hdulist[0], figure = plt.figure(0))
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
    ximgview(args.file_path)

