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


def xpselect(args):
    """Application for data subselection.

    We want to (loosely) model this on
    http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtselect.txt
    """
    evt_file_path = args.evfile
    assert(evt_file_path.endswith('.fits'))
    kwargs = args.__dict__
    xEventSelect(evt_file_path, **kwargs).select()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('evfile', type=str,
                        help='path to the input event file')
    parser.add_argument('--outfile', type=str, default=None,
                        help='path to the output binned FITS file')
    parser.add_argument('--ra', type=float, default=None,
                        help='RA of acceptance cone in decimal degrees')
    parser.add_argument('--dec', type=float, default=None,
                        help='Dec of acceptance cone in decimal degrees')
    parser.add_argument('--rad', type=float, default=None,
                        help='ROI radius in arcminutes')
    parser.add_argument('--tmin', type=float, default=None,
                        help='minimum time')
    parser.add_argument('--tmax', type=float, default=None,
                        help='maximum time')
    parser.add_argument('--phasemin', type=float, default=None,
                        help='minimum phase')
    parser.add_argument('--phasemax', type=float, default=None,
                        help='maximum phase')
    parser.add_argument('--emin', type=float, default=None,
                        help='minimum energy')
    parser.add_argument('--emax', type=float, default=None,
                        help='maximum energy')
    parser.add_argument('--phimin', type=float, default=None,
                        help='minimum azimuthal angle')
    parser.add_argument('--phimax', type=float, default=None,
                        help='maximum azimuthal angle')
    parser.add_argument('--mcsrcid', action='append', type=int, default=[],
                        help='the Monte Carlo source ID to select')
    parser.add_argument('--mc', action='store_true', default=False,
                        help='use Monte Carlo information for the selection')
    args = parser.parse_args()
    startmsg()
    xpselect(args)
