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
import os

from ximpol import XIMPOL_SRCMODEL
from ximpol.core.rand import xUnivariateGenerator
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
from ximpol.srcmodel.spectrum import power_law, xSourceSpectrum, constant


def _full_path(file_name):
    """
    """
    return os.path.join(XIMPOL_SRCMODEL, 'ascii', file_name)

SPEC_FILE_PATH = _full_path('Crab_PhaseResolvedSpectrum.txt')
PDEG_FILE_PATH = _full_path('Crab_Polarization_Degree.txt')
PANG_FILE_PATH = _full_path('Crab_Polarization_Angle.txt')

_phi0, _phi1, _idx, _idx_err, _norm, _norm_err = numpy.loadtxt(SPEC_FILE_PATH,
                                                               unpack=True)
_phi = 0.5*(_phi0 + _phi1)
pl_normalization = xInterpolatedUnivariateSpline(_phi, _norm, k=3)
pl_index = xInterpolatedUnivariateSpline(_phi, _idx, k=2)

_phi, _pol_angle = numpy.loadtxt(PANG_FILE_PATH, unpack=True)
_pol_angle = numpy.deg2rad(_pol_angle)
fmt = dict(xname=r'$\phi$', yname='Polarization angle [rad]')
pol_angle = xInterpolatedUnivariateSpline(_phi, _pol_angle, k=1, **fmt)

_phi, _pol_degree = numpy.loadtxt(PDEG_FILE_PATH, unpack=True)
_pol_degree /= 100.
fmt = dict(xname=r'$\phi$', yname='Polarization degree')
pol_degree = xInterpolatedUnivariateSpline(_phi, _pol_degree, k=1, **fmt)

def polarization_degree(E, t):
    return pol_degree(t)

def polarization_angle(E, t):
    return pol_angle(t)


ROI_MODEL = xROIModel(83.633083, 22.014500)

crab_ephemeris = xEphemeris(0., 29.8003951530036, -3.73414e-10, 1.18e-20)
crab_pulsar = xPeriodicPointSource('Crab pulsar', ROI_MODEL.ra, ROI_MODEL.dec,
                                   crab_ephemeris)
E = numpy.linspace(1, 10, 100)
t = numpy.linspace(0, 1, 100)
crab_pulsar.spectrum = power_law(pl_normalization, pl_index)
crab_pulsar.polarization_degree = polarization_degree
crab_pulsar.polarization_angle = polarization_angle

ROI_MODEL.add_source(crab_pulsar)


if __name__ == '__main__':
    print(ROI_MODEL)
    crab_pulsar.plot()
