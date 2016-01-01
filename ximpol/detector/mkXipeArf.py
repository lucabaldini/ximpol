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


import os
import numpy
from astropy.io import fits

import ximpol.__utils__
from ximpol import XIMPOL_DETECTOR, XIMPOL_IRF
from ximpol.__logging__ import logger
from ximpol.utils.xFunction1dTxtFile import xFunction1dTxtFile
from ximpol.fileio.xFitsDataFormatArf import xFitsDataFormatArf
from ximpol.detector.__XipeBaseline__ import *


def mkXipeArf():
    """ Create the .arf file for the XIPE baseline configuration.
    """
    outputFileName = 'xipe_%s.arf' % IRF_LABEL
    outputFilePath = os.path.join(XIMPOL_IRF, 'fits', outputFileName)
    if os.path.exists(outputFilePath):
        ximpol.__utils__.rm(outputFilePath)
    logger.info('Loading mirror effective area...')
    optsAeff = xFunction1dTxtFile(OPTS_AEFF_FILE_PATH, 'linear',
                                  xmin = ENERGY_MIN, xmax = ENERGY_MAX)
    logger.info('Loading detector quantum efficiency...')
    gpdEff = xFunction1dTxtFile(GPD_QEFF_FILE_PATH, 'linear',
                                xmin = ENERGY_MIN, xmax = ENERGY_MAX)
    logger.info('Evaluating effective area...')
    xipeAeff = optsAeff*gpdEff
    logger.info('Filling in arrays...')
    specresp = xipeAeff(ENERGY_MEAN)
    logger.info('Done, %d effectivea area values calculated.' % len(specresp))
    logger.info('Creating PRIMARY header and HDU...')
    primaryHeader = xFitsDataFormatArf.primaryHeader()
    print(repr(primaryHeader))
    primaryHdu = fits.PrimaryHDU(header = primaryHeader)
    logger.info('Creating SPECRESP header and HDU...')
    PRIMARY_HEADER_KWARGS['RESPFILE'] = outputFileName
    specrespHeader = xFitsDataFormatArf.specrespHeader(RESP_HEADER_COMMENTS,
                                                       **PRIMARY_HEADER_KWARGS)
    print(repr(specrespHeader))
    logger.info('Filling in SPECRESP data...')
    cols = xFitsDataFormatArf.specrespColumns([ENERGY_LO, ENERGY_HI, specresp])
    specrefHdu = fits.BinTableHDU.from_columns(cols, header = specrespHeader)
    logger.info('Writing output file...')
    hdulist = fits.HDUList([primaryHdu, specrefHdu])
    hdulist.info()
    hdulist.writeto(outputFilePath)
    logger.info('Done, bye!')
    return outputFilePath


if __name__ == '__main__':
    filePath = mkXipeArf()
    logger.info('Reading back file...')
    from ximpol.fileio.xInputArfFile import xInputArfFile
    f = xInputArfFile(filePath)
