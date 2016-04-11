#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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


import numpy
import sys
import os
import scipy
import scipy.signal

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger

#Not sure that I have the right ra and dec for the source, check!!
CYGX1_RA = 7.56
CYGX1_DEC = 6.59

MIN_ENERGY = 0.1
MAX_ENERGY = 15.

PL_INDEX = 2.3
#I intgrated the spline from 1 -10 keV
PL_NORM = 2.08

POL_DEGREE_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'CygX1_poldegree_model.txt')

POL_ANGLE_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'CygX1_polangle_model.txt')

FLUX_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'CygX1_flux_model.txt')


def polarization_degree(E, t, ra, dec):
    return pol_degree_spline(E)

def polarization_angle(E, t, ra, dec):
    return pol_angle_spline(E)

def energy_spectrum(E, t):
    return _energy_spectrum(E)


# Build the polarization degree as a function of the energy.
_energy, _pol_degree = numpy.loadtxt(POL_DEGREE_FILE_PATH, unpack=True)
_pol_degree /= 100.
print "Here are the energy and pol degree values"
print _energy
print
print _pol_degree
# Filter the data points to reduce the noise.
#_pol_degree = scipy.signal.wiener(_pol_degree, 5)
_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_pol_degree = _pol_degree[_mask]
fmt = dict(xname='Energy', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSpline(_energy, _pol_degree, k=1, **fmt)


# Build the polarization angle as a function of the energy.
_energy, _pol_angle = numpy.loadtxt(POL_ANGLE_FILE_PATH, unpack=True)
_pol_angle = numpy.deg2rad(_pol_angle)
print "Here are the energy and pol angle values"
print _energy
print
print _pol_angle
# Filter the data points to reduce the noise.
#_pol_angle = scipy.signal.wiener(_pol_angle, 2)
_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_pol_angle = _pol_angle[_mask]
fmt = dict(xname='Energy', yname='Polarization angle [rad]')
pol_angle_spline = xInterpolatedUnivariateSpline(_energy, _pol_angle, k=1, **fmt)

#"""
#put together the flux for the source starting from the txt file
_energy, _flux = numpy.loadtxt(FLUX_FILE_PATH, unpack=True)
_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
_energy = _energy[_mask]
_flux = _flux[_mask]
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
               yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')

_energy_spectrum = xInterpolatedUnivariateSpline(_energy, _flux, **fmt)
#"""

ROI_MODEL = xROIModel(CYGX1_RA, CYGX1_DEC)



src = xPointSource('Cyg-X1', ROI_MODEL.ra, ROI_MODEL.dec, energy_spectrum,
                   polarization_degree, polarization_angle)

ROI_MODEL.add_source(src)



if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    fig = plt.figure('Spectrum')
    _energy_spectrum.plot(show=False)
    fig = plt.figure('Polarization angle')
    pol_angle_spline.plot(show=False)
    fig = plt.figure('Polarization degree')
    pol_degree_spline.plot(show=False)
    plt.show()
