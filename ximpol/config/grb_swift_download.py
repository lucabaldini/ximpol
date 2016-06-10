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

import urllib2
import lxml.html
import os
import re

from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger

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

       min_obs_time: float
               minimum time required to the observation [seconds].
               Default: 22000 -> about 6 hours

       light_curve: str
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
            logger.info('Obs time < %f seconds!\nNot saving any data.'\
                        %min_obs_time)

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

def main():
    """
    """
    GRB_NAME = 'GRB 041223'
    flux_outfile = download_swift_grb_lc_file(GRB_NAME)
    rate_outfile = download_swift_grb_lc_file(GRB_NAME, light_curve='rate')
    
    grb_lc_ascii_file_list = []
    for grb in get_all_swift_grb_names():
        grb_lc_ascii_file = download_swift_grb_lc_file(grb)
        if type(grb_lc_ascii_file) is str:
            grb_lc_ascii_file_list.append(grb_lc_ascii_file)
    logger.info('Saved %i files.'%len(grb_lc_ascii_file_list))
    

if __name__=='__main__':
    main()
