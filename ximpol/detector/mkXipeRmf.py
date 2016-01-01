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
import scipy
from scipy import stats
#from stats import norm

import ximpol.__utils__
from ximpol import XIMPOL_DETECTOR, XIMPOL_IRF
from ximpol.__logging__ import logger
from ximpol.utils.xFunction1dTxtFile import xFunction1dTxtFile
from ximpol.utils.xFunction1d import xFunction1d
from ximpol.fileio.xFitsDataFormatRmf import xFitsDataFormatRmf
from ximpol.detector.__XipeBaseline__ import *
import matplotlib.pyplot as plt


def mkXipeRmf():
    """Create the .rmf file for the XIPE baseline configuration.
    """
    outputFileName = 'xipe_%s.rmf' % IRF_LABEL
    outputFilePath = os.path.join(XIMPOL_IRF, 'fits', outputFileName)
    if os.path.exists(outputFilePath):
        ximpol.__utils__.rm(outputFilePath)

    logger.info('Loading energy response function...')
    gpderes = xFunction1dTxtFile(GPD_ERES_FILE_PATH, 'linear')

    #Note that the eres is given only for energies between 2 and 8 keV whereas we computed the arf and mdf for 1-10 keV.

    logger.info('Filling in arrays...')



    gpderesresp = gpderes(ENERGY_MEAN)
    logger.info('Done, %d energy response values calculated.' %\
                len(gpderesresp))


    full_matrix_list = []
    channel_list = []
    Emin = []
    Emax = []

    x_range_min = 0.0
    x_range_max = 15.0

    for channel,energy in enumerate(ENERGY_MEAN):
        fwhm = gpderes(energy)
        rms = fwhm/2.3548


        #x_range = 4*rms
        #x_range_min = energy - x_range
        #x_range_max = energy + x_range

        #Filling in the list for EBOUNDS hhdulist, EMIN ad EMAX
        #accounding to what is listed here:ftp://legacy.gsfc.nasa.gov/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.ps page 15

        Emin.append( x_range_min)
        Emax.append( x_range_max)

        channel_list.append(channel)

        rv = stats.norm(loc=energy, scale=rms)
        x = numpy.linspace(x_range_min,x_range_max, NUM_CHANNELS)
        pdfsample = rv.pdf(x)
        matrix_sub = numpy.array([x,pdfsample])

        full_matrix_list.append(matrix_sub)

    F_CHAN = numpy.empty(len(ENERGY_MEAN))
    F_CHAN.fill(1)
    N_GRP = numpy.array([1])
    N_CHAN = numpy.empty(NUM_CHANNELS)
    N_CHAN.fill(NUM_CHANNELS)
    CHANNEL = numpy.array(channel_list)
    matrix = numpy.array(full_matrix_list)
    EMIN = numpy.array(Emin)
    EMAX = numpy.array(Emax)

    logger.info('Creating PRIMARY header and HDU...')
    primaryHeader = xFitsDataFormatRmf.primaryHeader()
    print(repr(primaryHeader))
    primaryHdu = fits.PrimaryHDU(header = primaryHeader)
    logger.info('Creating ERESRESP header and HDU...')
    PRIMARY_HEADER_KWARGS['ERESFILE'] = outputFileName
    eresrespHeader = xFitsDataFormatRmf.eresrespHeader(RESP_HEADER_COMMENTS,
                                                       **PRIMARY_HEADER_KWARGS)
    print(repr(eresrespHeader))
    logger.info('Filling in ERESRESP data...')

    cols = xFitsDataFormatRmf.eresrespColumns([ENERGY_LO, ENERGY_HI,
                                               N_GRP, F_CHAN,N_CHAN, matrix])
    eresrefHdu = fits.BinTableHDU.from_columns(cols, header = eresrespHeader)


    #Now I need to create the EBOUNDS header and col
    ereseboundsHeader = xFitsDataFormatRmf.ereseboundsHeader(RESP_HEADER_COMMENTS,
                                                       **PRIMARY_HEADER_KWARGS)
    print(repr(eresrespHeader))
    logger.info('Filling in ERESEBOUNS data...')

    eboundscols = xFitsDataFormatRmf.ereseboundsColumns([CHANNEL, EMIN,
                                                         EMAX])
    ereseboundsrefHdu = fits.BinTableHDU.from_columns(eboundscols, header = ereseboundsHeader)

    logger.info('Writing output file...')
    hdulist = fits.HDUList([primaryHdu, eresrefHdu, ereseboundsrefHdu])
    hdulist.info()
    hdulist.writeto(outputFilePath)
    logger.info('Done, bye!')
    return outputFilePath



if __name__ == '__main__':
    filePath = mkXipeRmf()
    logger.info('Reading back file...')
    from ximpol.fileio.xInputRmfFile import xInputRmfFile
    f = xInputRmfFile(filePath)
