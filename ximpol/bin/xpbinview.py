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


__description__ = 'Generic display interface to ximpol data products'


from astropy.io import fits

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg, abort
from ximpol.evt.binning import xBinnedCountSpectrum
from ximpol.evt.binning import xBinnedMap
from ximpol.evt.binning import xBinnedModulationCube
from ximpol.evt.binning import xBinnedLightCurve
from ximpol.evt.binning import xBinnedPhasogram


VIEW_DICT = {
    'PHA1' : xBinnedCountSpectrum,
    'LC'   : xBinnedLightCurve,
    'PHASG': xBinnedPhasogram,
    'CMAP' : xBinnedMap,
    'MCUBE': xBinnedModulationCube
}

def xpbinview(file_path):
    """Quick FITS image viewer.
    """
    try:
        binalg = fits.open(file_path)[0].header['BINALG']
    except Exception as e:
        abort('Could not determine file type (%s)' % e)
    VIEW_DICT[binalg](file_path).plot(show=True)


if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the input FITS image file')
    args = parser.parse_args()
    startmsg()
    xpbinview(args.file_path)
