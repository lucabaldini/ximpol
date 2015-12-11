#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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
from ximpol.__package__ import XIMPOL_DETECTOR, XIMPOL_IRF
from ximpol.__logging__ import logger
from ximpol.utils.xFunction1dTxtFile import xFunction1dTxtFile
from ximpol.fileio.xFitsDataFormatMrf import xFitsDataFormatMrf
from ximpol.detector.__XipeBaseline__ import *



def mkXipeMrf():
    """ Create the .mrf file for the XIPE baseline configuration.
    """
    outputFileName = 'xipe_%s.mrf' % IRF_LABEL
    outputFilePath = os.path.join(XIMPOL_IRF, 'fits', outputFileName)
    if os.path.exists(outputFilePath):
        ximpol.__utils__.rm(outputFilePath)
        
    logger.info('Loading modulation function...')
    gpdModf = xFunction1dTxtFile(GPD_MODF_FILE_PATH, 'linear')
    

    logger.info('Filling in arrays...')

    gpdmodfresp = gpdModf(ENERGY_MEAN)
    logger.info('Done, %d Modulation factor values calculated.' %\
                len(gpdmodfresp))
    logger.info('Creating PRIMARY header and HDU...')
    primaryHeader = xFitsDataFormatMrf.primaryHeader()
    print(repr(primaryHeader))
    primaryHdu = fits.PrimaryHDU(header = primaryHeader)
    logger.info('Creating MODFRESP header and HDU...')
    PRIMARY_HEADER_KWARGS['MODFFILE'] = outputFileName
    modfrespHeader = xFitsDataFormatMrf.modfrespHeader(RESP_HEADER_COMMENTS,
                                                       **PRIMARY_HEADER_KWARGS)
    print(repr(modfrespHeader))
    logger.info('Filling in MODFRESP data...')
    cols = xFitsDataFormatMrf.modfrespColumns([ENERGY_LO, ENERGY_HI,
                                               gpdmodfresp])
    modfrefHdu = fits.BinTableHDU.from_columns(cols, header = modfrespHeader)
    logger.info('Writing output file...')
    hdulist = fits.HDUList([primaryHdu, modfrefHdu])
    hdulist.info()
    hdulist.writeto(outputFilePath)
    logger.info('Done, bye!')
    return outputFilePath
    
    

if __name__ == '__main__':
    filePath = mkXipeMrf()
    logger.info('Reading back file...')
    from ximpol.fileio.xInputMrfFile import xInputMrfFile
    f = xInputMrfFile(filePath)
