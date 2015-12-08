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

from ximpol.__package__ import XIMPOL_DETECTOR



""" Basic configuration file with the XIPE characteristics, as listed in the
proposal to ESA.
"""

BASE_FOLDER = os.path.join(XIMPOL_DETECTOR, 'data')

OPTS_AEFF_FILE_NAME = 'aeff_optics_xipe_m4_x3.asc'
OPTS_AEFF_FILE_PATH = os.path.join(BASE_FOLDER, OPTS_AEFF_FILE_NAME)

GPD_QEFF_FILE_NAME = 'eff_hedme8020_1atm_1cm_cuts80p_be50um_p_x.asc'
GPD_QEFF_FILE_PATH = os.path.join(BASE_FOLDER, GPD_QEFF_FILE_NAME)

GPD_ERES_FILE_NAME = 'eres_fwhm_hedme8020_1atm_1cm.asc'
GPD_ERES_FILE_PATH = os.path.join(BASE_FOLDER, GPD_ERES_FILE_NAME)

GPD_MODF_FILE_NAME = 'modfact_hedme8020_1atm_1cm_mng.asc'
GPD_MODF_FILE_PATH = os.path.join(BASE_FOLDER, GPD_MODF_FILE_NAME)

GAS_MIXTURE = 'Ne/DME 80/20'
GAS_PRESSURE = 1.      # Atm
ABS_GAP_THICKNESS = 1. # cm
QUAL_CUT_EFFICIENCY = 0.8
WINDOW_MATERIAL = 'Be'
WINDOW_THICKNESS = 50. # um

ENERGY_MIN = 1.0       # keV
ENERGY_MAX = 9.9       # keV
ENERGY_STEP = 0.01     # keV
ENERGY_LO = numpy.arange(ENERGY_MIN, ENERGY_MAX, ENERGY_STEP)
ENERGY_HI = numpy.append(ENERGY_LO[1:], ENERGY_MAX)
ENERGY_MEAN = 0.5*(ENERGY_LO + ENERGY_HI)


PRIMARY_HEADER_KWARGS = {
    'TELESCOP': 'XIPE',
    'INSTRUME': 'GPD',
    'NAXIS2'  : len(ENERGY_MEAN),
}

SPECRESP_HEADER_COMMENTS = [
    'Gas mixture: %s' % GAS_MIXTURE,
    'Pressure: %.3f Atm' % GAS_PRESSURE,
    'Absorption gap: %.3f cm' % ABS_GAP_THICKNESS,
    'Quality cut efficiency: %.3f' % QUAL_CUT_EFFICIENCY,
    'Window: %s, %d um' % (WINDOW_MATERIAL, WINDOW_THICKNESS)
]
