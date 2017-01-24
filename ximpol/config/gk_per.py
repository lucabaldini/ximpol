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

from ximpol import XIMPOL_CONFIG
from ximpol.srcmodel.roi import xPeriodicPointSource, xEphemeris, xROIModel
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.units_ import keV2erg


GK_PER_RA = 52.80005
GK_PER_DEC = 43.9043
GK_PER_PERIOD = 531.3


def _full_path(file_name):
    """Convenience function to retrieve the relevant files.
    """
    return os.path.join(XIMPOL_CONFIG, 'ascii', file_name)

"""Parse the input file with the energy spectrum.

This has been converted into a text file by Niccolo` starting from an XSPEC
model that Domitilla provided. There's also a simplified version of the model
in the same folder called `gk_per_spectrum_simple.txt` but since from the
standpoint of the simulation that doesn't make any difference, we're
using the right one.

Note that the normaliztion of the integral energy flux between 2 and 8 keV
is provided into a separate file (see below) so that we need to integrate
this to have the conversion factors.
"""
# Read in the file and create the spline.
spectrum_file_path = _full_path('gk_per_spectrum_full.txt')
_emin, _emax, _flux = numpy.loadtxt(spectrum_file_path, unpack=True)
_e = numpy.sqrt(_emin*_emax)
_fmt = dict(xname='Energy', xunits='keV', yname='Differential flux',
            yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
spectrum_spline = xInterpolatedUnivariateSplineLinear(_e, _flux, **_fmt)
# And we need a spline for the energy flux as well, so that we can integrate it
# and evaluate the average energy flux. This comes out around
# 1.0214032297e-10 erg/s/cm2, so it's in the right ballpark---good.
_fmt = dict(xname='Energy', xunits='keV', yname='Differential energy flux',
            yunits='cm$^{-2}$ s$^{-1}$')
espectrum_spline =  xInterpolatedUnivariateSplineLinear(_e, _e*_flux, **_fmt)
average_eflux = keV2erg(espectrum_spline.integral(2., 8.))


"""Parse the phasogram file.

This includes the integral energy flux between 2 and 8 keV, in erg cm^-2 s^-1
in 10 phase bins.

Despite the fact that the input values are most likely the average fluxes
within each phase bin, we just add 0.05 and consider them as the values
at the bin center. In addition, in order not to have discontinuities at the
boundaries, we  set the values at 0 and 1 as the average values of the first
and last point. We might want to be smarter in the long run, but this should
be a sensible first step.
"""
phasogram_file_path = _full_path('gk_per_flux_2_8_keV_ph10.txt')
_bin, _phase, _eflux, _dummy = numpy.loadtxt(phasogram_file_path, unpack=True)
# Duplicate the values at phi=0 for phi=1.
_phase = numpy.append(_phase, 1.)
_eflux = numpy.append(_eflux, _eflux[0])
_fmt = dict(xname='Phase', yname='Integral energy flux 2--8 keV',
            yunits='erg s$^{-1}$ cm$^{-2}$')
phasogram_spline = xInterpolatedUnivariateSpline(_phase, _eflux, k=3, **_fmt)


"""Parse the input polarization degree and angle.

This is provided in a text file in ten phase bins, and we exactly the same
magic that we do for the integral flux above.
"""
pol_degree_file_path = _full_path('gk_per_pol_degree_angle_ph10.txt')
_dummy, _phase, _degree, _angle = numpy.loadtxt(pol_degree_file_path,
                                                unpack=True)
#_degree /= 100.
# Duplicate the values at phi=0 for phi=1.
_phase = numpy.append(_phase, 1.)
_degree = numpy.append(_degree, _degree[0])
_angle = numpy.append(_angle, _angle[0])
_fmt = dict(xname='Phase', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSpline(_phase, _degree, k=3, **_fmt)
_fmt = dict(xname='Phase', yname='Polarization angle [$^\\circ$]')
pol_angle_spline = xInterpolatedUnivariateSpline(_phase, _angle, k=3, **_fmt)


def energy_spectrum(E, phase):
    """Return the energy spectrum.
    """
    return phasogram_spline(phase)/average_eflux*spectrum_spline(E)

def polarization_degree(E, phase, ra, dec):
    """Return the polarization degree.
    """
    return pol_degree_spline(phase)

def polarization_angle(E, phase, ra, dec):
    """Return the polarization angle.
    """
    return numpy.radians(pol_angle_spline(phase))


ephem = xEphemeris(0., 1./GK_PER_PERIOD)


"""We have all the ingredients, can define the ROI.

Note that, since the input XSPEC model is already including the interstellar
abrosorption, we are setting the column density explicitely to 0 here.
"""
ROI_MODEL = xROIModel(GK_PER_RA, GK_PER_DEC)
gk_per = xPeriodicPointSource('GK Per', GK_PER_RA, GK_PER_DEC, energy_spectrum,
                              polarization_degree, polarization_angle, ephem,
                              column_density=0., redshift=0.)
ROI_MODEL.add_source(gk_per)



if __name__ == '__main__':
    print(ROI_MODEL)
    from ximpol.utils.matplotlib_ import pyplot as plt
    plt.figure()
    spectrum_spline.plot(show=False, logx=True, logy=True)
    plt.axis([0.75, 11, 1e-10, 1e-2])
    plt.figure()
    phasogram_spline.plot(show=False)
    plt.figure()
    pol_degree_spline.plot(show=False)
    plt.figure()
    pol_angle_spline.plot(show=False)
    plt.show()
