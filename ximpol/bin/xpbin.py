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


__description__ = 'Bin the data'


import os
import numpy
from astropy.io import fits

from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.logging_ import logger, startmsg
from ximpol.irf.base import xPrimaryHDU, update_header
from ximpol.evt.binning import xColDefsSpectrum, SPECTRUM_HEADER_SPECS


def xpbin(file_path):
    """Application to bin the data.
    """
    assert(file_path.endswith('.fits'))
    logger.info('Opening input file %s...' % file_path)
    hdu_list = fits.open(file_path)
    hdu_list.info()
    evtdata = hdu_list['EVENTS'].data

    # Horrible---need to put the EBOUND information in the event file.
    _num_chans = 256
    # And also the GTI
    _elapsed_time = evtdata['TIME'][-1] - evtdata['TIME'][0]

    _binning = numpy.linspace(-0.5, _num_chans - 0.5, _num_chans)
    n, bins, patches = plt.hist(evtdata['PHA'], bins=_binning)

    output_file_path = file_path.replace('.fits', '.pha')
    primary_hdu = xPrimaryHDU()
    data = [numpy.arange(_num_chans),
            n/_elapsed_time,
            numpy.sqrt(n)/_elapsed_time
    ]
    cols = xColDefsSpectrum(data)
    spectrum_hdu = fits.BinTableHDU.from_columns(cols)
    update_header(spectrum_hdu, SPECTRUM_HEADER_SPECS)
    hdulist = fits.HDUList([primary_hdu, spectrum_hdu])
    hdulist.info()
    hdulist.writeto(output_file_path, clobber=True)
    logger.info('Event list written to %s...' % output_file_path)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the path to the input event file')
    args = parser.parse_args()
    startmsg()
    xpbin(args.file_path)
