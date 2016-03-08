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


from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.select import xEventSelect


"""Command-line switches.
"""
import argparse
PARSER = argparse.ArgumentParser(description=__description__)
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


def xpselect(file_path, **kwargs):
    """Application for data subselection.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtselect.txt
    """
    assert(file_path.endswith('.fits'))
    return xEventSelect(file_path, **kwargs).select()


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpselect(args.evfile, **args.__dict__)
