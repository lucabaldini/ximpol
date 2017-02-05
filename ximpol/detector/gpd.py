#!/usr/bin/env python
#
# Copyright (C) 2017, the ximpol team.
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

from ximpol import XIMPOL_DETECTOR
from ximpol.utils.logging_ import logger

def _full_path(file_name):
    return os.path.join(XIMPOL_DETECTOR, 'data', file_name)


BE_DENSITY = 1.85 # g cm^{-3}

N = 6.022140857e23
R = 0.08206 
H_MASS = 1.008
HE_MASS = 4.002602
C_MASS = 12.011
O_MASS = 15.999
DME_MASS = 2*C_MASS + 6*H_MASS + O_MASS


def perfect_gas_molar_volume(temperature, pressure):
    """Return the molar volume (in cm^3) at a given temperature and pressure.

    This is using the law of perfect gases.

    Arguments
    ---------
    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in atm.
    """
    return R*(temperature + 273.15)/pressure*1000.

def perfect_gas_density(mass, temperature, pressure):
    """Return the density (in g cm^{-3}) for a perfect gas of a given 
    molecular mass at a given temperature and pressure.

    Arguments
    ---------
    mass : float
        The molecular mass of the gas.

    temperature : float
        The temperature in degrees C.

    pressure : float
        The pressure in atm.
    """
    return mass/perfect_gas_molar_volume(temperature, pressure)

def window_transparency(thickness=0.0050, separator=','):
    """Calculate the transparency of a Be window and write the output to file.

    Arguments
    ---------
    thickness : float
        The window thickness (in cm).

    separator : string
        The column separator for the output file.
    """
    data = numpy.loadtxt(_full_path('xcom_be_xsec.txt'), unpack=True)
    energy, coher_xsec, incoher_xsec, photo_xsec = data
    # Convert MeV to keV
    energy *= 1000.
    lambda_inv = photo_xsec*BE_DENSITY
    transparence = numpy.exp(-lambda_inv*thickness)
    file_path = _full_path('window_transparency_be_%dum.csv' %\
                           (thickness*10000))
    logger.info('Writing output file %s...' % file_path)
    output_file = open(file_path, 'w')
    output_file.write('#Energy%sTransparency\n' % separator)
    output_file.write('#[keV]%s[]\n' % separator)
    for _e, _trans in zip(energy, transparence):
        output_file.write('%.2f%s%.6f\n' % (_e, separator, _trans))
    output_file.close()
    logger.info('Done.')

def quantum_efficiency(thickness=1.0, he_pp=0.2, dme_pp=0.8, T=0.,
                       cut_efficiency=1., separator=','):
    """Calculate the quantum efficiency of a Ne/DME gass cell.
    """
    he_data = numpy.loadtxt(_full_path('xcom_he_xsec.txt'), unpack=True)
    energy, he_coher_xsec, he_incoher_xsec, he_photo_xsec = he_data
    dme_data = numpy.loadtxt(_full_path('xcom_dme_xsec.txt'), unpack=True)
    energy, dme_coher_xsec, dme_incoher_xsec, dme_photo_xsec = dme_data
    # Convert MeV to keV
    energy *= 1000.
    lambda_inv = he_photo_xsec*perfect_gas_density(HE_MASS, T, he_pp) +\
                 dme_photo_xsec*perfect_gas_density(DME_MASS, T, dme_pp)
    efficiency = 1 - numpy.exp(-lambda_inv*thickness)
    efficiency *= cut_efficiency
    file_path = _full_path('gpd_quantum_efficiency_he%d_dme%d_%dmm.csv' %\
                           (he_pp*100, dme_pp*100, thickness*10))
    if cut_efficiency < 1.:
        file_path = file_path.replace('.csv', '_cut%d.csv' %\
                                      (cut_efficiency*100))
    logger.info('Writing output file %s...' % file_path)
    output_file = open(file_path, 'w')
    output_file.write('#Energy%sEfficiency\n' % separator)
    output_file.write('#[keV]%s[]\n' % separator)
    for _e, _eff in zip(energy, efficiency):
        output_file.write('%.4f%s%.5f\n' % (_e, separator, _eff))
    output_file.close()
    logger.info('Done.')

    

if __name__ == '__main__':
    window_transparency()
    quantum_efficiency()
    quantum_efficiency(cut_efficiency=0.8)

