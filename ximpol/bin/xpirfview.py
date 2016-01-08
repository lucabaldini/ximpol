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


__description__ = 'Quick IRF viewer'


from ximpol.utils.matplotlib_ import pyplot as plt
from astropy.io import fits

from ximpol.utils.logging_ import startmsg, abort
from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.mrf import xModulationFactor
from ximpol.irf.psf import xPointSpreadFunction
from ximpol.irf.rmf import xEnergyDispersion


CLASS_DICT = {
    'arf': xEffectiveArea,
    'mrf': xModulationFactor,
    'rmf': xEnergyDispersion,
    'psf': xPointSpreadFunction
}


def xpirfview(file_path):
    """Quick FITS image viewer.
    """
    file_ext = file_path.split('.').pop()
    if not file_ext in CLASS_DICT.keys():
        abort('Unrecognized file extension (.%s)' % file_ext)
    irf = CLASS_DICT[file_ext](file_path)
    irf.plot()


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the input FITS IRF file')
    args = parser.parse_args()
    startmsg()
    xpirfview(args.file_path)
