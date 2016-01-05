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

from scipy import stats
from astropy.io import fits

from ximpol import XIMPOL_DETECTOR, XIMPOL_IRF
from ximpol.utils.logging_ import logger
from ximpol.utils.os_ import rm
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.irf.base import xPrimaryHDU
from ximpol.irf.arf import SPECRESP_HEADER_SPECS, xColDefsSPECRESP
from ximpol.irf.mrf import MODFRESP_HEADER_SPECS, xColDefsMODFRESP
from ximpol.irf.rmf import EBOUNDS_HEADER_SPECS, MATRIX_HEADER_SPECS,\
    xColDefsMATRIX, xColDefsEBOUNDS



"""Basic configuration file with the XIPE characteristics, as listed in the
proposal to ESA.
"""

IRF_LABEL = 'baseline'
BASE_FOLDER = os.path.join(XIMPOL_DETECTOR, 'data')
OPT_AEFF_FILE_PATH = os.path.join(BASE_FOLDER,
                    'aeff_optics_xipe_m4_x3.asc')
GPD_QEFF_FILE_PATH = os.path.join(BASE_FOLDER,
                    'eff_hedme8020_1atm_1cm_cuts80p_be50um_p_x.asc')
GPD_ERES_FILE_PATH = os.path.join(BASE_FOLDER,
                    'eres_fwhm_hedme8020_1atm_1cm.asc')
GPD_MODF_FILE_PATH = os.path.join(BASE_FOLDER,
                    'modfact_hedme8020_1atm_1cm_mng.asc')
GAS_MIXTURE = 'Ne/DME 80/20'#
GAS_PRESSURE = 1.           # Atm
ABS_GAP_THICKNESS = 1.      # cm
QUAL_CUT_EFFICIENCY = 0.8   #
WINDOW_MATERIAL = 'Be'      #
WINDOW_THICKNESS = 50.      # um
ENERGY_MIN = 1.0            # keV
ENERGY_MAX = 10.0           # keV
ENERGY_STEP = 0.01          # keV
NUM_CHANNELS = 1024         #
E_CHAN_OFFSET = 0.          # keV
E_CHAN_SLOPE = 0.01         # keV/channel
ENERGY_LO = numpy.arange(ENERGY_MIN, ENERGY_MAX, ENERGY_STEP)
ENERGY_HI = numpy.append(ENERGY_LO[1:], ENERGY_MAX)
ENERGY_MEAN = 0.5*(ENERGY_LO + ENERGY_HI)


XIPE_HEADER_SPEC = [
    ('TELESCOP', 'XIPE' , 'mission/satellite name'),
    ('INSTRUME', 'GPD'  , 'instrument/detector name'),
    ('DETNAM'  , 'NONE' , 'specific detector name in use')
]

SPECRESP_HEADER_SPECS += XIPE_HEADER_SPEC
MODFRESP_HEADER_SPECS += XIPE_HEADER_SPEC
EBOUNDS_HEADER_SPECS += XIPE_HEADER_SPEC

RESP_HEADER_COMMENTS = [
    'Gas mixture: %s' % GAS_MIXTURE,
    'Pressure: %.3f Atm' % GAS_PRESSURE,
    'Absorption gap: %.3f cm' % ABS_GAP_THICKNESS,
    'Quality cut efficiency: %.3f' % QUAL_CUT_EFFICIENCY,
    'Window: %s, %d um' % (WINDOW_MATERIAL, WINDOW_THICKNESS)
]


def make_arf():
    """Write the XIPE effective area response function.
    """
    logger.info('Creating XIPE effective area fits file...')
    output_file_name = 'xipe_%s.arf' % IRF_LABEL
    output_file_path = os.path.join(XIMPOL_IRF, 'fits', output_file_name)
    if os.path.exists(output_file_path):
        rm(output_file_path)
    logger.info('Loading mirror effective area from %s...' % OPT_AEFF_FILE_PATH)
    _x, _y = numpy.loadtxt(OPT_AEFF_FILE_PATH, unpack=True)
    opt_aeff = xInterpolatedUnivariateSplineLinear(_x, _y, ENERGY_MIN,
                                                   ENERGY_MAX)
    logger.info('Loading quantum efficiency from %s...' % GPD_QEFF_FILE_PATH)
    _x, _y = numpy.loadtxt(GPD_QEFF_FILE_PATH, unpack=True)
    gpd_eff = xInterpolatedUnivariateSplineLinear(_x, _y, ENERGY_MIN,
                                                  ENERGY_MAX)
    logger.info('Evaluating effective area...')
    aeff = opt_aeff*gpd_eff
    logger.info('Filling in arrays...')
    specresp = aeff(ENERGY_MEAN)
    logger.info('Creating PRIMARY HDU...')
    primary_hdu = xPrimaryHDU()
    print(repr(primary_hdu.header))
    logger.info('Creating SPECRESP HDU...')
    cols = xColDefsSPECRESP([ENERGY_LO, ENERGY_HI, specresp])
    specresp_hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, value, comment in SPECRESP_HEADER_SPECS:
        specresp_hdu.header.set(keyword, value, comment)
    for comment in RESP_HEADER_COMMENTS:
        specresp_hdu.header['COMMENT'] = comment
    print(repr(specresp_hdu.header))
    logger.info('Writing output file %s...' % output_file_path)
    hdulist = fits.HDUList([primary_hdu, specresp_hdu])
    hdulist.info()
    hdulist.writeto(output_file_path)
    logger.info('Done.')


def make_mrf():
    """
    """
    logger.info('Creating XIPE effective area fits file...')
    output_file_name = 'xipe_%s.mrf' % IRF_LABEL
    output_file_path = os.path.join(XIMPOL_IRF, 'fits', output_file_name)
    if os.path.exists(output_file_path):
        rm(output_file_path)
    logger.info('Loading modulation factor from %s...' % GPD_MODF_FILE_PATH)
    _x, _y = numpy.loadtxt(GPD_MODF_FILE_PATH, unpack=True)
    modf = xInterpolatedUnivariateSplineLinear(_x, _y, ENERGY_MIN, ENERGY_MAX)
    logger.info('Filling in arrays...')
    modfresp = modf(ENERGY_MEAN)
    logger.info('Creating PRIMARY HDU...')
    primary_hdu = xPrimaryHDU()
    print(repr(primary_hdu.header))
    logger.info('Creating MODFRESP HDU...')
    cols = xColDefsMODFRESP([ENERGY_LO, ENERGY_HI, modfresp])
    modfresp_hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, value, comment in MODFRESP_HEADER_SPECS:
        modfresp_hdu.header.set(keyword, value, comment)
    for comment in RESP_HEADER_COMMENTS:
        modfresp_hdu.header['COMMENT'] = comment
    print(repr(modfresp_hdu.header))
    logger.info('Writing output file %s...' % output_file_path)
    hdulist = fits.HDUList([primary_hdu, modfresp_hdu])
    hdulist.info()
    hdulist.writeto(output_file_path)
    logger.info('Done.')


def make_rmf():
    """

    ftp://legacy.gsfc.nasa.gov/caldb/docs/memos/cal_gen_92_002/cal_gen_92_002.ps
    page 15
    """
    output_file_name = 'xipe_%s.rmf' % IRF_LABEL
    output_file_path = os.path.join(XIMPOL_IRF, 'fits', output_file_name)
    if os.path.exists(output_file_path):
        rm(output_file_path)
    logger.info('Loading energy dispersion from %s...' % GPD_ERES_FILE_PATH)
    _x, _y = numpy.loadtxt(GPD_ERES_FILE_PATH, unpack=True)
    edisp_fwhm = xInterpolatedUnivariateSplineLinear(_x, _y, ENERGY_MIN,
                                                     ENERGY_MAX)
    logger.info('Creating PRIMARY HDU...')
    primary_hdu = xPrimaryHDU()
    print(repr(primary_hdu.header))
    logger.info('Creating MATRIX HDU...')
    nrows = len(ENERGY_LO)
    ngrp = numpy.ones(nrows)
    fchan = numpy.zeros(nrows)
    nchan = numpy.array([NUM_CHANNELS]*nrows, 'i')
    matrix = numpy.zeros((0, NUM_CHANNELS))
    ch = numpy.arange(NUM_CHANNELS)
    for energy, rms in zip(ENERGY_MEAN, edisp_fwhm(ENERGY_MEAN)/2.358):
        mean_chan = int((energy - E_CHAN_OFFSET)/E_CHAN_SLOPE)
        rms_chan = rms/E_CHAN_SLOPE
        rv = stats.norm(loc=mean_chan, scale=rms_chan)
        matrix = numpy.vstack([matrix, rv.pdf(ch)])
    cols = xColDefsMATRIX([ENERGY_LO, ENERGY_HI, ngrp, fchan, nchan, matrix])
    matrix_hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, value, comment in MATRIX_HEADER_SPECS:
        matrix_hdu.header.set(keyword, value, comment)
    print(repr(matrix_hdu.header))
    logger.info('Creating EBOUNDS HDU...')
    ch = numpy.arange(NUM_CHANNELS)
    emin = ch*E_CHAN_SLOPE + E_CHAN_OFFSET
    emax = (ch + 1)*E_CHAN_SLOPE + E_CHAN_OFFSET
    cols = xColDefsEBOUNDS([ch, emin, emax])
    ebounds_hdu = fits.BinTableHDU.from_columns(cols)
    for keyword, value, comment in EBOUNDS_HEADER_SPECS:
        ebounds_hdu.header.set(keyword, value, comment)
    print(repr(ebounds_hdu.header))
    logger.info('Writing output file %s...' % output_file_path)
    hdulist = fits.HDUList([primary_hdu, matrix_hdu, ebounds_hdu])
    hdulist.info()
    hdulist.writeto(output_file_path)
    logger.info('Done.')



if __name__ == '__main__':
    make_arf()
    make_mrf()
    make_rmf()
