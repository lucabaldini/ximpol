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

import lxml.html
import urllib2
import numpy
import os
import re

from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear

PATH_0 = 'http://www.swift.ac.uk'
PATH_1 = os.path.join(PATH_0,'xrt_curves')
ALL_LC_PATH = os.path.join(PATH_1,'allcurvesflux.php')
GRB_ALL_NAME = re.compile('[a-zA-Z]+\s[a-zA-Z]*\s*\d+[A-Z]*')

curve = {'flux':'flux.qdp',
         'rate':'curve.qdp',
         'details_rate':'curve2.qdp'}

def get_all_swift_grb_names():
    """Returns the list of all the Swift GRBs extracted from
       http://www.swift.ac.uk/xrt_curves/allcurvesflux.php
    """
    grb_list = []
    f = urllib2.urlopen(ALL_LC_PATH)
    for line in f:
        if "class='grb'>" in line:
            m = GRB_ALL_NAME.search(line)
            grb_list.append(m.group())
    return grb_list

def download_swift_grb_lc_file(grb_name, min_obs_time=22000, \
                               light_curve='flux'):
    """Downloads an ascii file with the light curve of a given GRB 
       observed by Swift.
       
       Arguments                                                               
       ---------                                                               
       grb_name: str
               GRB name, as given (including spaces) in 
               http://www.swift.ac.uk/xrt_curves/allcurvesflux.php

       min_obs_time: float, optional
               minimum time required to the observation [seconds].
               Default: 22000 -> about 6 hours

       light_curve: str, optional
               Defines the kind of curve: if flux(t) or rate(t).
               Default: 'flux' -> [erg cm^{-2} s^{-1}]
               Other options: 'rate' -> [s^{-1}]
    """
    outpath = os.path.join(XIMPOL_CONFIG,'ascii','grb_swift')
    if not os.path.exists(outpath):
        os.makedirs(outpath)
    outfile = os.path.join(outpath,'grb_%s_%s_Swift.dat'\
                           %(grb_name.replace(' ','').replace('.','-'),
                             light_curve))
    if os.path.exists(outfile):
        logger.info('Already saved %s'%outfile)
        return outfile
    else:
        lc_url, lc = get_lc_url(grb_name,light_curve)
        data = lc.read().split()
        last_time = float(data[len(data)-6])
        if check_obs_time(last_time,min_obs_time):
            ff = open(outfile,'w')
            lc = urllib2.urlopen(lc_url)
            ff.write(lc.read())
            ff.close()
            logger.info('Saving %s'%outfile)
            return outfile
        elif not check_obs_time(last_time,min_obs_time):
            logger.info('Obs time for %s < %f seconds!\nNot saving any data.'\
                        %(grb_name,min_obs_time))

def check_obs_time(last_time, min_obs_time):
    """check if the time of the observation is greater than a referece time

       Arguments                                                               
       ---------  
       min_obs_time: float
               minimum time required to the observation [seconds].

       last_time: str
               The time of the last event after the trigger.
    """
    if last_time >= min_obs_time:
        return True
    else: 
        return False

def get_lc_url(grb_name, light_curve):
    """get the url of the light curve ascii file for a given GRB

       Arguments                                                               
       ---------  
       grb_name: str
               GRB name, as given (including spaces) in 
               http://www.swift.ac.uk/xrt_curves/allcurvesflux.php
       light_curve: str
               Defines the kind of curve: if flux(t) or rate(t).
    """
    f = urllib2.urlopen(ALL_LC_PATH)
    curve_file = curve[light_curve]
    lc_url, lc = '', ''
    for line in f:
        if grb_name in line:
            dom =  lxml.html.fromstring(line)
            for link in dom.xpath('//a/@href'):
                lc_url = PATH_0+link+curve_file
                lc = urllib2.urlopen(lc_url)
    return lc_url, lc    

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

def parse_light_curve(file_path, num_bins=100, num_min_data=5):
    """Read the light-curve data points and make sense out of them.

    Here we make a weighted average of the input data over a logarithmic
    time binning, then we interpolate in log space to fill the gaps in the
    input data, and finally we create a spline in linear space.
    """
    t, tp, tm, f, fp, fm = parse_light_curve_data(file_path)
    if len(t) >= num_min_data:
        binning = numpy.logspace(numpy.log10(t[0]), numpy.log10(t[-1]), \
                                 num_bins)
        tave = [t[0]]
        fave = [f[0]]
        for _tmin, _tmax in zip(t[:-1], t[1:]):
            #_tave = (_tmin*_tmax)**0.5
            _tave = _tmin
            _mask = (t >= _tmin)*(t <= _tmax)
            if numpy.count_nonzero(_mask) > 1:
                _weights = numpy.power(0.5*(fp + fm)[_mask], -2.)
                _fave = numpy.average(f[_mask], weights=_weights)
                _fave = f[_mask][0]
                tave.append(_tave)
                fave.append(_fave)
        """
        for _tmin, _tmax in zip(binning[:-1], binning[1:]):
            _tave = (_tmin*_tmax)**0.5
            _mask = (t >= _tmin)*(t <= _tmax)
            if numpy.count_nonzero(_mask) > 1:
                _weights = numpy.power(0.5*(fp + fm)[_mask], -2.)
                _fave = numpy.average(f[_mask], weights=_weights)
                tave.append(_tave)
                fave.append(_fave)"""
        tave = numpy.log10(numpy.array(tave))
        fave = numpy.log10(numpy.array(fave))
        spline = xInterpolatedUnivariateSplineLinear(tave, fave)
        #tave = numpy.linspace(tave[0], tave[-1], num_bins)
        tave = numpy.linspace(tave[0], tave[-1], len(t))
        fave = spline(tave)
        tave = numpy.power(10., tave)
        fave = numpy.power(10., fave)
        fmt = dict(xname='Time', xunits='s',
                   yname='Energy integral flux 0.3-10 keV',
                   yunits='erg cm$^{-2}$ s$^{-1}$')
        return xInterpolatedUnivariateSplineLinear(tave, fave, **fmt)
    else:
        logger.info('Data points < %i !\nNot producing any light curve.'\
                        %num_min_data)
        return 0

def main():
    """
    """
    grb_lc_ascii_file_list = []
    for grb in get_all_swift_grb_names():
        grb_lc_ascii_file = download_swift_grb_lc_file(grb)
        if type(grb_lc_ascii_file) is str:
            grb_lc_ascii_file_list.append(grb_lc_ascii_file)
    logger.info('Saved %i files.'%len(grb_lc_ascii_file_list))
    

if __name__=='__main__':
    main()
