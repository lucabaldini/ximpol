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

"""References:
Slowikowska A. et al., MNRAS 397, 103-23

Valid range (MJD)       : 52944--52975
Epoch, t0 (MJD)         : 52960.000000296
nu0 (Hz)                : 29.8003951530036
nudot(10^-10 Hz s^-1)   : -3.73414
nudddot (10^-20 Hz s^-2): 1.18
"""


import numpy
import scipy.signal
from scipy.optimize import curve_fit
import os

from ximpol import XIMPOL_CONFIG
from ximpol.core.rand import xUnivariateGenerator
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
from ximpol.srcmodel.spectrum import power_law


def _full_path(file_name):
    """Convenience function to retrieve the relevant files.
    """
    return os.path.join(XIMPOL_CONFIG, 'ascii', file_name)

# Grab all the relevant files.
SPEC_FILE_PATH = _full_path('Crab_PhaseResolvedSpectrum.txt')
PDEG_FILE_PATH = _full_path('Crab_Polarization_Degree.txt')
PANG_FILE_PATH = _full_path('Crab_Polarization_Angle.txt')

# Read the file with the phase-resolved spectral parameters.
_phi0, _phi1, _index, _index_err, _norm,\
    _norm_err = numpy.loadtxt(SPEC_FILE_PATH, unpack=True)
_phi = 0.5*(_phi0 + _phi1)

# Build the PL normalization as a function of the phase.
fmt = dict(xname='Pulsar phase', yname='PL normalization',
           yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
pl_normalization_spline = xInterpolatedUnivariateSpline(_phi, _norm, k=3, **fmt)

# Fit the PL index as a function of the phase with a sinusoid.
def ffit(x, *pars):
    """Fit function for the PL index as function of the phase.
    """
    return pars[0] + pars[1]*numpy.cos(pars[2]*x + pars[3])

p0 = [1., 1., 1., 1.]
popt, pcov = curve_fit(ffit, _phi, _index, p0, _index_err)

# And use tho best-fit parameters to build a proper spline.
fmt = dict(xname='Pulsar phase', yname='PL index')
pl_index_spline = xInterpolatedUnivariateSpline(_phi, ffit(_phi, *popt), k=2,
                                                **fmt)

# Build the actual energy spectrum.
energy_spectrum = power_law(pl_normalization_spline, pl_index_spline)

# Build the polarization angle as a function of the phase.
_phi, _pol_angle = numpy.loadtxt(PANG_FILE_PATH, unpack=True)
_pol_angle = numpy.deg2rad(_pol_angle)
# Filter the data points to reduce the noise.
_pol_angle = scipy.signal.wiener(_pol_angle, 7)
fmt = dict(xname='Pulsar phase', yname='Polarization angle [rad]')
pol_angle_spline = xInterpolatedUnivariateSpline(_phi, _pol_angle, k=1, **fmt)

# Mind that you have to wrap this into a function to be used.
# This could be done in a more general way with some sort of library function.
def polarization_degree(E, t, ra, dec):
    return pol_degree_spline(t)

# Build the polarization degree as a function of the phase.
_phi, _pol_degree = numpy.loadtxt(PDEG_FILE_PATH, unpack=True)
_pol_degree /= 100.
# Filter the data points to reduce the noise.
_pol_degree = scipy.signal.wiener(_pol_degree, 5)
fmt = dict(xname='Pulsar phase', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSpline(_phi, _pol_degree, k=1, **fmt)

# And, again, this needs to be wrapped into a function.
def polarization_angle(E, t, ra, dec):
    return pol_angle_spline(t)


ROI_MODEL = xROIModel(83.633083, 22.014500)
crab_ephemeris = xEphemeris(0., 29.8003951530036, -3.73414e-10, 1.18e-20)
crab_pulsar = xPeriodicPointSource('Crab pulsar', ROI_MODEL.ra, ROI_MODEL.dec,
                                   energy_spectrum, polarization_degree,
                                   polarization_angle, crab_ephemeris)
ROI_MODEL.add_source(crab_pulsar)


if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    plt.figure()
    pl_index_spline.plot(show=False)
    plt.figure()
    pl_normalization_spline.plot(show=False)
    plt.figure()
    pol_angle_spline.plot(show=False)
    plt.figure()
    pol_degree_spline.plot(show=False)
    plt.show()
