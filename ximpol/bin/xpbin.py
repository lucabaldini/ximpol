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
BIN_ALGORITHMS = BIN_ALG_DICT.keys()
BIN_ALGORITHMS.sort()


def xpbin(args):
    """Application to bin the data.
    """
    evt_file_path = args.evfile
    algorithm = args.algorithm
    assert(evt_file_path.endswith('.fits'))
    kwargs = args.__dict__
    kwargs.pop('evfile')
    kwargs.pop('algorithm')
    BIN_ALG_DICT[algorithm](evt_file_path, **kwargs).bin_()


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('--evfile', type=str,
                        help='the path to the input event file')
    parser.add_argument('--algorithm', choices=BIN_ALGORITHMS, required=True,
                        help='the binning mode')
    parser.add_argument('--outfile', type=str, default=None,
                        help='the output binned FITS file')
    args = parser.parse_args()
    startmsg()
    xpbin(args)
