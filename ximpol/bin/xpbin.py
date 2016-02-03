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


__description__ = 'Bin event data in different flavors'


import os
import numpy
from astropy.io import fits

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.binning import xEventBinningPHA1
from ximpol.evt.binning import xEventBinningLC
from ximpol.evt.binning import xEventBinningCMAP


BIN_ALG_DICT = {'PHA1': xEventBinningPHA1,
                'LC'  : xEventBinningLC,
                'CMAP': xEventBinningCMAP
}
BIN_ALGS = BIN_ALG_DICT.keys()
BIN_ALGS.sort()
TBIN_ALGS = ['FILE', 'LIN', 'LOG', 'SNR']
PRJCTS = ['AIT', 'ZEA', 'ARC', 'CAR', 'GLS', 'MER', 'NCP', 'SIN', 'STG', 'TAN']
COORD_SYS = ['CEL', 'GAL']


def xpbin(args):
    """Application to bin the data.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtbin.txt
    """
    evt_file_path = args.evfile
    algorithm = args.algorithm
    assert(evt_file_path.endswith('.fits'))
    kwargs = args.__dict__
    kwargs.pop('evfile')
    kwargs.pop('algorithm')
    print kwargs
    BIN_ALG_DICT[algorithm](evt_file_path, **kwargs).bin_()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('--evfile', type=str,
                        help='path to the input event file')
    parser.add_argument('--algorithm', choices=BIN_ALGS, required=True,
                        help='the binning algorithm')
    parser.add_argument('--outfile', type=str, default=None,
                        help='path to the output binned FITS file')
    parser.add_argument('--tbinalg', choices=TBIN_ALGS, default='LIN',
                        help='time binning specification')
    parser.add_argument('--tstart', type=float, default=None,
                        help='start time for LIN/LOG time binning')
    parser.add_argument('--tstop', type=float, default=None,
                        help='stop time for LIN/LOG time binning')
    parser.add_argument('--tbins', type=int, default=100,
                        help='number of bins for LIN/LOG time binning')
    parser.add_argument('--tbinfile', type=str, default=None,
                        help='path to the optional time bin definition file')
    parser.add_argument('--nxpix', type=int, default=256,
                        help='number of horizontal pixels in the output image')
    parser.add_argument('--nypix', type=int, default=256,
                        help='number of vertical pizels in the output image')
    parser.add_argument('--binsz', type=float, default=1.,
                        help='the pixel size of the output image in arcseconds')
    parser.add_argument('--xref', type=float, default=None,
                        help='the horizontal position of the image center')
    parser.add_argument('--yref', type=float, default=None,
                        help='the vertical position of the image center')
    parser.add_argument('--proj', choices=PRJCTS, default='TAN',
                        help='coordinate projection')
    args = parser.parse_args()
    startmsg()
    xpbin(args)
