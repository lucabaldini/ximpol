#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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
import os

from ximpol.srcmodel.gabs import xpeInterstellarAbsorptionModel
from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.utils.units_ import keV2erg

RA = 203.973
DEC = -34.2960

SPIN = 1          #spin of the black hole (0 or 1)
INCLINATION = 60  #degree (fixed)
INT_POL_DEG = 2   #intrinsic polarization degree (0, 2 or 4)
INT_POL_ANG = 0   #intrinsic polarization angle (0 or 90)
FLUX = 1e-11      #erg/s/cm2
GAL_NH = 0.14e+21 #cm-2

if INT_POL_ANG is 90:
    H = 10        #height in gravitational radii (3 or 10)
else:
    H = 3

file_name = 'lamp_pol_%d_%d_0%02d_0000_%d_%02d.dat' %(SPIN, INCLINATION, H,
            INT_POL_DEG, INT_POL_ANG)
file_path = os.path.join(XIMPOL_CONFIG, 'ascii', file_name)

def parse(file_path, emin=1., emax=10.):
    """Parse the input file with the git complete spectral and polarization model.
    """
    logger.info('Parsing input file %s...' % file_path)
    energy, flux, degree, angle = numpy.loadtxt(file_path, unpack=True,
                                                        usecols = (0,1,5,6))
    if INT_POL_ANG is 0:
        angle += 90.
    _mask = (energy >= emin)*(energy <= emax)
    return energy[_mask], flux[_mask], degree[_mask], angle[_mask]
    #return energy[_mask], flux[_mask], numpy.full(degree[_mask].shape, INT_POL_DEG/100.), numpy.full(angle[_mask].shape, INT_POL_ANG+90.)

def scale_flux(energy, flux, emin=2., emax=8.):
    ism_model = xpeInterstellarAbsorptionModel()
    ism_trans = ism_model.transmission_factor(GAL_NH)
    flux *= ism_trans(energy)
    norm = xInterpolatedUnivariateSplineLinear(energy,energy*flux).integral(
                                                                emin, emax)
    norm = keV2erg(norm)
    scale_factor =  FLUX/norm
    return flux*scale_factor

def energy_spectrum(E, t):
    return spectral_model(E)

def polarization_degree(E, t, ra, dec):
    return pol_degree(E)

def polarization_angle(E, t, ra, dec):
    return numpy.radians(pol_angle(E))

energy, flux, degree, angle = parse(file_path)
flux = scale_flux(energy, flux)
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
           yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
spectral_model = xInterpolatedUnivariateSplineLinear(energy, flux, **fmt)
fmt = dict(xname='Energy', xunits='keV', yname='Polarization degree')
pol_degree = xInterpolatedUnivariateSplineLinear(energy, degree, **fmt)
fmt = dict(xname='Energy', xunits='keV', yname='Polarization angle',
           yunits='$^\\circ$')
pol_angle = xInterpolatedUnivariateSplineLinear(energy, angle, **fmt)

ROI_MODEL = xROIModel(RA, DEC)

source = xPointSource('NGC 1365', RA, DEC, energy_spectrum,
                      polarization_degree, polarization_angle)

ROI_MODEL.add_source(source)

def display():
    """Display the source model.
    """
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    fig = plt.figure('Energy spectrum')
    spectral_model.plot(show=False, logy=True)
    fig = plt.figure('Polarization degree')
    pol_degree.plot(show=False)
    fig = plt.figure('Polarization angle')
    pol_angle.plot(show=False)
    plt.show()

if __name__=='__main__':
    display()
