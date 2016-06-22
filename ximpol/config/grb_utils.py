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
import os

from ximpol.utils.logging_ import logger
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


def get_grb_spec_index(data_file):
    f = open(data_file)
    for line in f:
        if 'INDEX = ' in line:
            l = line.split()
            index = l[len(l)-1]
    return float(index)

def get_grb_position(data_file):
    f = open(data_file)
    for line in f:
        if 'RA = ' in line:
            l = line.split()
            ra = l[len(l)-1]
        if 'DEC = ' in line:
            l = line.split()
            dec = l[len(l)-1]
    return float(ra), float(dec)

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

def parse_light_curve_original(file_path, num_bins=150):
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
        else:
            _t1 = numpy.log10(t[t <= _tmin][-1])
            _t2 = numpy.log10(t[t >= _tmax][0])
            _f1 = numpy.log10(f[t <= _tmin][-1])
            _f2 = numpy.log10(f[t >= _tmax][0])
            _t = numpy.log10(_tave)
            _f = _f1 + (_f2 - _f1)*(_t - _t1)/(_t2 - _t1)
            _fave = 10**(_f)
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

def parse_light_curve(file_path, num_bins=100, num_min_data=5):
    """Read the light-curve data points and make sense out of them.
    """
    t, tp, tm, f, fp, fm = parse_light_curve_data(file_path)
    return xInterpolatedUnivariateSplineLinear(t, f)
    print t
    if len(t) >= num_min_data:
        tave = [t[0]]
        fave = [f[0]]
        for _tmin, _tmax in zip(t[:-1], t[1:]):
            _tave = _tmin
            _mask = (t >= _tmin)*(t <= _tmax)
            if numpy.count_nonzero(_mask) > 1:
                _fave = f[_mask][0]
                tave.append(_tave)
                fave.append(_fave)
        print tave
        tave = numpy.log10(numpy.array(tave))
        fave = numpy.log10(numpy.array(fave))
        spline = xInterpolatedUnivariateSplineLinear(tave, fave)
        tave = numpy.linspace(tave[0], tave[-1], len(t))
        fave = spline(tave)
        tave = numpy.power(10., tave)
        fave = numpy.power(10., fave)
        fmt = dict(xname='Time', xunits='s',
                   yname='Integral energy flux 0.3--10 keV',
                   yunits='erg cm$^{-2}$ s$^{-1}$')
        return xInterpolatedUnivariateSplineLinear(tave, fave, **fmt)
    else:
        logger.info('Data points < %i !\nNot producing any light curve.'\
                        %num_min_data)
        return None

def main():
    """Test the script Retriving RA, Dec and Index for GRB 130427A 
    """
    from ximpol.utils.matplotlib_ import pyplot as plt
    grb_name = 'GRB 050219A'
    from ximpol.config.grb_swift_download import download_swift_grb_lc_file
    file_path = download_swift_grb_lc_file(grb_name)
    grb_ra, grb_dec = get_grb_position(file_path)
    grb_index = get_grb_spec_index(file_path)
    logger.info('Retriving information for %s:'%grb_name)
    logger.info('\tPosition: RA = %fdeg, Dec = %fdeg'%(grb_ra,grb_dec))
    logger.info('\tSpectral Index = %f '%grb_index)
    grb_lc = parse_light_curve_original(file_path)
    t, tp, tm, f, fp, fm = parse_light_curve_data(file_path)
    plt.errorbar(t, f, yerr=[fm, fp], fmt='o')
    grb_lc.plot(logx=True, logy=True)
    


if __name__=='__main__':
    main()
