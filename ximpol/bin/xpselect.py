#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


__description__ = 'Data selection interface'


import os

from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.subselect import xEventSelect


"""Command-line switches.
"""
import argparse
import ast

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('evfile', type=str,
                    help='path to the input event file')
PARSER.add_argument('--outfile', type=str, default=None,
                    help='path to the output binned FITS file')
PARSER.add_argument('--ra', type=float, default=None,
                    help='RA of acceptance cone in decimal degrees')
PARSER.add_argument('--dec', type=float, default=None,
                    help='Dec of acceptance cone in decimal degrees')
PARSER.add_argument('--rad', type=float, default=None,
                    help='ROI radius in arcminutes')
PARSER.add_argument('--tmin', type=float, default=None,
                    help='minimum time')
PARSER.add_argument('--tmax', type=float, default=None,
                    help='maximum time')
PARSER.add_argument('--phasemin', type=float, default=None,
                    help='minimum phase')
PARSER.add_argument('--phasemax', type=float, default=None,
                    help='maximum phase')
PARSER.add_argument('--emin', type=float, default=None,
                    help='minimum energy')
PARSER.add_argument('--emax', type=float, default=None,
                    help='maximum energy')
PARSER.add_argument('--phimin', type=float, default=None,
                    help='minimum azimuthal angle')
PARSER.add_argument('--phimax', type=float, default=None,
                    help='maximum azimuthal angle')
PARSER.add_argument('--mcsrcid', action='append', type=int, default=[],
                    help='the Monte Carlo source ID to select')
PARSER.add_argument('--mc', action='store_true', default=False,
                    help='use Monte Carlo information for the selection')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')


def xpselect(file_path, **kwargs):
    """Application for data subselection.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtselect.txt
    """
    assert(file_path.endswith('.fits'))
    event_select = xEventSelect(file_path, **kwargs)
    outfile = event_select.get('outfile')
    if os.path.exists(outfile) and not event_select.get('clobber'):
        logger.info('Output file %s already exists.' % outfile)
        logger.info('Remove the file or set "clobber = True" to overwite it.')
    else:
        event_select.select()
    return outfile


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpselect(args.evfile, **args.__dict__)
