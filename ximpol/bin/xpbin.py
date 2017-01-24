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

from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.binning import xEventBinningPHA1
from ximpol.evt.binning import xEventBinningLC
from ximpol.evt.binning import xEventBinningPHASG
from ximpol.evt.binning import xEventBinningCMAP
from ximpol.evt.binning import xEventBinningMCUBE


BIN_ALG_DICT = {
    'PHA1' : xEventBinningPHA1,
    'LC'   : xEventBinningLC,
    'PHASG': xEventBinningPHASG,
    'CMAP' : xEventBinningCMAP,
    'MCUBE': xEventBinningMCUBE
}
BIN_ALGS = BIN_ALG_DICT.keys()
BIN_ALGS.sort()
TBIN_ALGS = ['FILE', 'LIN', 'LOG']
EBIN_ALGS = ['FILE', 'LIN', 'LOG', 'EQP', 'LIST']
PRJCTS = ['AIT', 'ZEA', 'ARC', 'CAR', 'GLS', 'MER', 'NCP', 'SIN', 'STG', 'TAN']
COORD_SYS = ['CEL', 'GAL']


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('evfile', type=str,
                    help='path to the input event file')
PARSER.add_argument('--algorithm', choices=BIN_ALGS, required=True,
                    help='the binning algorithm')
PARSER.add_argument('--outfile', type=str, default=None,
                    help='path to the output binned FITS file')
PARSER.add_argument('--tbinalg', choices=TBIN_ALGS, default='LIN',
                    help='time binning specification')
PARSER.add_argument('--tstart', type=float, default=None,
                    help='start time for LIN/LOG time binning')
PARSER.add_argument('--tstop', type=float, default=None,
                    help='stop time for LIN/LOG time binning')
PARSER.add_argument('--tbins', type=int, default=100,
                    help='number of bins for LIN/LOG time binning')
PARSER.add_argument('--tbinfile', type=str, default=None,
                    help='path to the optional time bin definition file')
PARSER.add_argument('--phasebins', type=int, default=50,
                    help='number of bins for phase binning')
PARSER.add_argument('--nxpix', type=int, default=256,
                    help='number of horizontal pixels in the output image')
PARSER.add_argument('--nypix', type=int, default=256,
                    help='number of vertical pizels in the output image')
PARSER.add_argument('--binsz', type=float, default=2.5,
                    help='the pixel size of the output image in arcseconds')
PARSER.add_argument('--xref', type=float, default=None,
                    help='the horizontal position of the image center')
PARSER.add_argument('--yref', type=float, default=None,
                    help='the vertical position of the image center')
PARSER.add_argument('--proj', choices=PRJCTS, default='TAN',
                    help='coordinate projection')
PARSER.add_argument('--ebinalg', choices=EBIN_ALGS, default='LIN',
                    help='energy binning specification')
PARSER.add_argument('--emin', type=float, default=1.,
                    help='minimum energy for LIN/LOG energy binning')
PARSER.add_argument('--emax', type=float, default=10.,
                    help='maximum energy for LIN/LOG energy binning')
PARSER.add_argument('--ebins', type=int, default=5,
                    help='number of bins for LIN/LOG energy binning')
PARSER.add_argument('--ebinfile', type=str, default=None,
                    help='path to the optional energy bin definition file')
PARSER.add_argument('--ebinning', type=ast.literal_eval, default=None,
                    help='the list containing the bin edges')
PARSER.add_argument('--phibins', type=int, default=75,
                    help='number of bins for LIN/LOG phi binning')
PARSER.add_argument('--mcsrcid', type=ast.literal_eval, default=None,
                    help='the list of MC source ID to be considered as source')
PARSER.add_argument('--mc', type=ast.literal_eval, choices=[True, False],
                    default=False,
                    help='use Monte Carlo information for binning')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')


def xpbin(file_path, **kwargs):
    """Application to bin the data.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtbin.txt
    """
    assert(file_path.endswith('.fits'))
    event_binning = BIN_ALG_DICT[kwargs['algorithm']](file_path, **kwargs)
    outfile = event_binning.get('outfile')
    if os.path.exists(outfile) and not event_binning.get('clobber'):
        logger.info('Output file %s already exists.' % outfile)
        logger.info('Remove the file or set "clobber = True" to overwite it.')
    else:
        event_binning.bin_()
    return outfile


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpbin(args.evfile, **args.__dict__)
