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
from ximpol.fileio.xFitsDataFormatArf import xFitsDataFormatArf



def mkXipeArf(emin = 1., emax = 9.9, estep = 0.01):
    """ 
    """
    outputFileName = 'xipe_proposal.arf'
    outputFilePath = os.path.join(XIMPOL_IRF, 'fits', outputFileName)
    if os.path.exists(outputFilePath):
        ximpol.__utils__.rm(outputFilePath)
    
    baseFolder = os.path.join(XIMPOL_DETECTOR, 'data')
    logger.info('Loading mirror effective area...')
    optsAeffFileName = 'aeff_optics_xipe_m4_x3.asc'
    optsAeffFilePath = os.path.join(baseFolder, optsAeffFileName)
    optsAeff = xFunction1dTxtFile(optsAeffFilePath, 'linear',
                                  xmin = emin, xmax = emax)

    logger.info('Loading detector quantum efficiency...')
    gpdEffFileName = 'eff_hedme8020_1atm_1cm_cuts80p_be50um_p_x.asc'
    gpdEffFilePath = os.path.join(baseFolder, gpdEffFileName)
    gpdEff = xFunction1dTxtFile(gpdEffFilePath, 'linear',
                                  xmin = emin, xmax = emax)

    logger.info('Evaluating effective area...')
    xipeAeff = optsAeff*gpdEff

    logger.info('Filling in arrays...')
    elo = numpy.arange(emin, emax, estep)
    ehi = numpy.append(elo[1:], emax)
    emean = 0.5*(elo + ehi)
    specresp = xipeAeff(emean)
    nbins = len(specresp)
    logger.info('Done, %d effectivea area values calculated.' % nbins)

    logger.info('Creating PRIMARY header and HDU...')
    primaryHeader = xFitsDataFormatArf.primaryHeader()
    print(repr(primaryHeader))
    primaryHdu = fits.PrimaryHDU(header = primaryHeader)
    
    logger.info('Creating SPECRESP header and HDU...')
    kwspecs = {
        'TELESCOP': 'XIPE',
        'INSTRUME': 'GPD',
        'NAXIS2'  : nbins,
        'RESPFILE': outputFileName
    }
    comments = [
        'Gas mixture: Ne/DME 80/20',
        'Pressure: 1 Atm',
        'Absorption gap: 1 cm',
        'Quality cut efficiency: 80%',
        'Window: Be, 50 um'
    ]
    specrespHeader = xFitsDataFormatArf.specrespHeader(comments, **kwspecs)
    print(repr(specrespHeader))

    logger.info('Filling in SPECRESP data...')
    cols = xFitsDataFormatArf.specrespColumns([elo, ehi, specresp])
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
