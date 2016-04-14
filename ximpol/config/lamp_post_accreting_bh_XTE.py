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

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger


"""In the so called 'lamp-post' model, the unpolarized emission from a primary
source illuminates the accretion disk of a black hole (BH), where it is
reprocessed. Part of it is emitted toward the observer, with a polarization
degree and angle depending on the exact geometry of the model.

The main parameters of the 'lamp post' emission model are:
- the height of the primary source above to the accretion disk, in units of the
gravitational radius
- the observer inclination angle
- the angular momentum per unit mass of the black hole in units of the
gravitational radius. Spin = 0 correspond to a static BH (Schwarzchild metric),
Spin = 1 to an extremely rotating one (Kerr metric).

Here we are using spectrum, polarization degree and polarization angle
for a model of the transient black hole system XTE J1650-500 presented in
Dovciak et al. 'Light bending scenario for accreting black holes in x-ray
polarimetry' (2011), provided by Fabio Muleri.

The required info are stored as columns in several .dat files. Each file
corresponds to a different choice of the main input parameters.
"""


data_file_name = 'lamp_pol_1_30_004_reduced.dat'
data_file_path = os.path.join (XIMPOL_CONFIG, 'ascii', data_file_name)


def parse(file_path, emin=1., emax=15., flux_scale=4.3576):
    """Parse the input file with the complete spectral and polarization model.
    """
    logger.info('Parsing input file %s...' % file_path)
    energy, flux, angle, degree = numpy.loadtxt(file_path, unpack=True)
    flux *= flux_scale
    _mask = (energy >= emin)*(energy <= emax)
    return energy[_mask], flux[_mask], angle[_mask], degree[_mask]


energy, flux, degree, angle = parse(data_file_path)
fmt = dict(xname='Energy', xunits='keV', yname='Flux',
           yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
spectral_model = xInterpolatedUnivariateSplineLinear(energy, flux, **fmt)
fmt = dict(xname='Energy', xunits='keV', yname='Polarization degree')
pol_degree = xInterpolatedUnivariateSplineLinear(energy, degree, **fmt)
fmt = dict(xname='Energy', xunits='keV', yname='Polarization angle',
           yunits='$^\\circ$')
pol_angle = xInterpolatedUnivariateSplineLinear(energy, angle, **fmt)


def energy_spectrum(E, t):
    return spectral_model(E)

def polarization_degree(E, t, ra, dec):
    return pol_degree(E)

def polarization_angle(E, t, ra, dec):
    return numpy.radians(pol_angle(E))

RA = 336.7183
DEC = -03.4271
ROI_MODEL = xROIModel(RA, DEC)

source = xPointSource('J1650-500', RA, DEC, energy_spectrum,
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
