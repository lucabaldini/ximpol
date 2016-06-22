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
from ximpol.srcmodel.spectrum import power_law, int_eflux2pl_norm
from ximpol.srcmodel.polarization import constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG


GRB_RA = 173.1362
GRB_DEC = 27.7129
MIN_ENERGY = 0.3
MAX_ENERGY = 10.
PL_INDEX = 1.66


def parse_light_curve_data(file_path):
    """Parse the ASCII file with the GRB light curve data from Swift.
    """
    t = []
    tp = []
    tm = []
    f = []
    fp = []
    fm = []
    for line in open(file_path):
        try:
            _t, _tp, _tm, _f, _fp, _fm = [float(item) for item in line.split()]
            t.append(_t)
            tp.append(_tp)
            tm.append(_tm)
            f.append(_f)
            fp.append(_fp)
            fm.append(-_fm)
        except:
            pass
    t = numpy.array(t)
    tp = numpy.array(tp)
    tm = numpy.array(tm)
    f = numpy.array(f)
    fp = numpy.array(fp)
    fm = numpy.array(fm)
    return t, tp, tm, f, fp, fm

def parse_light_curve(file_path, num_bins=150):
    """Read the light-curve data points and make sense out of them.

    Here we make a weighted average of the input data over a logarithmic
    time binning, then we interpolate in log space to fill the gaps in the
    input data, and finally we create a spline in linear space.
    """
    t, tp, tm, f, fp, fm = parse_light_curve_data(file_path)
    binning = numpy.logspace(numpy.log10(t[0]), numpy.log10(t[-1]), num_bins)
    tave = [t[0]]
    fave = [f[0]]
    for _tmin, _tmax in zip(binning[:-1], binning[1:]):
        _tave = (_tmin*_tmax)**0.5
        _mask = (t >= _tmin)*(t <= _tmax)
        if numpy.count_nonzero(_mask) > 1:
            _weights = numpy.power(0.5*(fp + fm)[_mask], -2.)
            _fave = numpy.average(f[_mask], weights=_weights)
            tave.append(_tave)
            fave.append(_fave)
    tave = numpy.log10(numpy.array(tave))
    fave = numpy.log10(numpy.array(fave))
    spline = xInterpolatedUnivariateSplineLinear(tave, fave)
    tave = numpy.linspace(tave[0], tave[-1], num_bins)
    fave = spline(tave)
    tave = numpy.power(10., tave)
    fave = numpy.power(10., fave)
    fmt = dict(xname='Time', xunits='s',
               yname='Energy integral flux 0.3-10 keV',
               yunits='erg cm$^{-2}$ s$^{-1}$')
    return xInterpolatedUnivariateSplineLinear(tave, fave, **fmt)


ROI_MODEL = xROIModel(GRB_RA, GRB_DEC)

lc_file_path = os.path.join(XIMPOL_CONFIG, 'ascii/GRB130427_Swift.dat')
integral_flux_spline = parse_light_curve(lc_file_path)
scale_factor = int_eflux2pl_norm(1, MIN_ENERGY, MAX_ENERGY, PL_INDEX)
fmt = dict(yname='PL normalization', yunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')
pl_normalization_spline = integral_flux_spline.scale(scale_factor, **fmt)

"""For the polarization degree we are literally making up something
that decreases with time.
"""
_t = integral_flux_spline.x
_p = 0.6*(1. - (_t/integral_flux_spline.xmax())**0.1)
fmt = dict(xname='Time', xunits='s', yname='Polarization degree')
pol_degree_spline = xInterpolatedUnivariateSplineLinear(_t, _p, **fmt)

def energy_spectrum(E, t):
    return pl_normalization_spline(t)*numpy.power(E, -PL_INDEX)

def polarization_degree(E, t, ra, dec):
    return pol_degree_spline(t)

polarization_angle = constant(numpy.radians(28.))

grb = xPointSource('GRB', GRB_RA, GRB_DEC, energy_spectrum, polarization_degree,
                   polarization_angle,
                   min_validity_time=integral_flux_spline.xmin(),
                   max_validity_time=integral_flux_spline.xmax())

def sampling_time(tstart, tstop):
    return numpy.logspace(numpy.log10(tstart), numpy.log10(tstop), 100)

grb.sampling_time = sampling_time


ROI_MODEL.add_source(grb)


if __name__=='__main__':
    from ximpol.utils.matplotlib_ import pyplot as plt
    from ximpol.utils.matplotlib_ import save_current_figure
    from ximpol import XIMPOL_DOC
    output_folder = os.path.join(XIMPOL_DOC, 'figures', 'showcase')

    fig = plt.figure()
    integral_flux_spline.plot(logx=True, logy=True, show=False)
    if True:
        save_current_figure('grb130427_swift_input_lc', output_folder, False)
    fig = plt.figure()
    pl_normalization_spline.plot(logx=True, logy=True, show=False)
    fig = plt.figure()
    pol_degree_spline.plot(logx=True, logy=False, show=False)
    plt.show()
