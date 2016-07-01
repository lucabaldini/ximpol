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


import os
import numpy
from scipy.interpolate import RectBivariateSpline

from ximpol import XIMPOL_CONFIG
from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedBivariateSpline
from ximpol.utils.units_ import keV2erg, erg2keV
from ximpol.utils.matplotlib_ import pyplot as plt


RA = 0.
DEC = 0.
PERIOD = 11.005
COLUMN_DENSITY = 0.#1.36e22
NUM_PHASE_BINS = 9
NUM_ENERGY_BINS = 35
MAXIMUM_ENERGY = 15.


file_path = os.path.join(XIMPOL_CONFIG, 'ascii', 'J1708_90_60_05_034.dat')
input_file = open(file_path)

def read_line():
    return [float(item) for item in input_file.next().strip('\n').split()]

chi, xi = read_line()
print chi, xi



flux_data = numpy.zeros((NUM_PHASE_BINS, NUM_ENERGY_BINS))
pol_degree_data = numpy.zeros((NUM_PHASE_BINS, NUM_ENERGY_BINS))
pol_angle_data = numpy.zeros((NUM_PHASE_BINS, NUM_ENERGY_BINS))
phase_min_data = []
phase_max_data = []
energy_data = []

for i in range(NUM_PHASE_BINS):
    index, phase_min, phase_max = read_line()
    phase_min_data.append(phase_min)
    phase_max_data.append(phase_max)
    for j in range(NUM_ENERGY_BINS):
        energy_min, energy_max, flux, pol_degree, pol_angle = read_line()
        energy_min = 10.**energy_min
        energy_max = 10.**energy_max
        energy_mean = numpy.sqrt(energy_min*energy_max)
        pol_angle = numpy.radians(pol_angle)
        flux_data[i, j] = flux
        pol_degree_data[i, j] = pol_degree
        pol_angle_data[i, j] = pol_angle
        if i == 0:
            energy_data.append(energy_mean)

phase_min_data = numpy.array(phase_min_data)
phase_max_data = numpy.array(phase_max_data)
phase_mean_data = 0.5*(phase_min_data + phase_max_data)/(2*numpy.pi)
energy_data = numpy.array(energy_data)

_x = phase_mean_data
_y = energy_data[energy_data<MAXIMUM_ENERGY]
_z = flux_data[:,energy_data<MAXIMUM_ENERGY]
_fmt = dict(xname='Phase', xunits='rad', yname='Energy', yunits='keV',
            zname='Differential flux', zunits='cm$^-2$ s$^{-1}$ keV$^{-1}$')
energy_spectrum_spline = xInterpolatedBivariateSpline(_x, _y, _z, kx=3, ky=1,
                                                      **_fmt)

_z = pol_degree_data[:,energy_data<MAXIMUM_ENERGY]
_fmt = dict(xname='Phase', xunits='rad', yname='Energy', yunits='keV',
            zname='Polarization degree')
polarization_degree_spline = xInterpolatedBivariateSpline(_x, _y, _z,
                                                          kx=3, ky=1, **_fmt)

_z = pol_angle_data[:,energy_data<MAXIMUM_ENERGY]
_fmt = dict(xname='Phase', xunits='rad', yname='Energy', yunits='keV',
            zname='Polarization angle', zunits='rad')
polarization_angle_spline = xInterpolatedBivariateSpline(_x, _y, _z, kx=3, ky=1,
                                                         **_fmt)

    
def energy_spectrum(E, phase):
    """
    """
    return energy_spectrum_spline.__call__(phase, E)


def polarization_degree(E, phase, ra, dec):
    """
    """
    return polarization_degree_spline.__call__(phase, E)


def polarization_angle(E, phase, ra, dec):
    """
    """
    return polarization_angle_spline.__call__(phase, E)


ephem = xEphemeris(0., 1./PERIOD)


"""We have all the ingredients, can define the ROI.
"""
ROI_MODEL = xROIModel(RA, DEC)
axp_1708 = xPeriodicPointSource('J1708', RA, DEC, energy_spectrum,
                                polarization_degree, polarization_angle, ephem,
                                column_density=COLUMN_DENSITY, redshift=0.)
ROI_MODEL.add_source(axp_1708)


    

if __name__ == '__main__':
    plt.figure()
    energy_spectrum_spline.plot(show=False)

    plt.figure()
    polarization_degree_spline.plot(show=False)

    plt.figure()
    polarization_angle_spline.plot(show=False)
    
    plt.figure()
    _x = numpy.linspace(0, 1, 100)
    _y = energy_spectrum(2., _x)
    plt.plot(_x, _y)

    plt.figure()
    _x = numpy.linspace(1, 10, 100)
    for ph in numpy.linspace(0, 0.1, 2):
        _y = polarization_degree(_x, ph, 0, 0)
        plt.plot(_x, _y, label='Phase = %.2f' % ph)
    plt.legend(loc='upper right')

    plt.figure()
    _x = numpy.linspace(1, 10, 100)
    for ph in numpy.linspace(0, 0.1, 2):
        _y = polarization_angle(_x, ph, 0, 0)
        plt.plot(_x, _y, label='Phase = %.2f' % ph)
    plt.legend(loc='upper right')

    plt.figure()
    _x = numpy.linspace(0, 1, 100)
    for e in numpy.linspace(1, 10, 10):
        _y = polarization_degree(e, _x, 0, 0)
        plt.plot(_x, _y, label='Energy = %.2f keV' % e)
    plt.legend(loc='upper right')
        
    plt.show()
