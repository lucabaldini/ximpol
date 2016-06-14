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

PATH_0 = 'http://www.swift.ac.uk'
PATH_CURV = os.path.join(PATH_0,'xrt_curves')
PATH_POS = os.path.join(PATH_0,'xrt_positions')
PATH_UNENH_POS = os.path.join(PATH_0,'xrt_unenh_positions')
PATH_SPEC = os.path.join(PATH_0,'xrt_spectra')
ALL_LC_PATH = os.path.join(PATH_CURV,'allcurvesflux.php')
GRB_ALL_NAME = re.compile('[a-zA-Z]+\s[a-zA-Z]*\s*\d+[A-Z]*')
PHOTON_INDEX = re.compile('\>\d\.\d+\s')
RA_DEC = re.compile('\([-]*\d+\.\d+')
#NH = re.compile('\d+\.*\d*')

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
    grb_list.remove('GRB 120711B')
    return grb_list

def download_swift_grb_lc_file(grb_name, min_obs_time=21600, \
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
        logger.info('Already saved file for %s'%grb_name)
        return outfile
    else:
        lc_url, lc = get_lc_url(grb_name,light_curve)
        data = lc.read().split()
        last_time = float(data[len(data)-6])
        if check_obs_time(last_time,min_obs_time):
            index = get_grb_spec_index_from_site(grb_name)
            ra, dec = get_grb_position_from_site(grb_name)
            ff = open(outfile,'w')
            lc = urllib2.urlopen(lc_url)
            ff.write('%s INDEX = %f\n'%(grb_name,index))
            ff.write('%s RA = %f\n'%(grb_name,ra))
            ff.write('%s DEC = %f\n'%(grb_name,dec))
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

def get_trigger_number(grb_name):
    trigger = ''
    f = urllib2.urlopen(ALL_LC_PATH)
    for line in f:
        if grb_name in line:
            dom = lxml.html.fromstring(line)
            for link in dom.xpath('//a/@href'):
                trigger = link.replace('xrt_curves','')
                trigger = trigger.replace('/','')
            break
    return trigger

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
    curve_file = curve[light_curve]
    lc_url, lc = '', ''
    trigger = get_trigger_number(grb_name)
    lc_url = os.path.join(PATH_CURV,trigger,curve_file)
    lc = urllib2.urlopen(lc_url)
    return lc_url, lc

def get_grb_spec_index_from_site(grb_name):
    """Extract the photon index from Swift web site for a given GRB.
       If more than one photon index, it returns the late time spectrum index.

       Arguments
       ---------
       grb_name: str
               GRB name, as given (including spaces) in
               http://www.swift.ac.uk/xrt_curves/allcurvesflux.php
    """
    trigger = get_trigger_number(grb_name)
    spec_url = os.path.join(PATH_SPEC,trigger)
    spec_page = urllib2.urlopen(spec_url)
    grb_spec_index = 0.
    for line in spec_page:
        if 'Photon index' in line:
            m = PHOTON_INDEX.search(line)
            if m == None:
                grb_spec_index = 2.
                logger.info('Spectrum not accurate. Default index set to 2.')
            else:
                grb_spec_index = float(m.group().replace('>',''))
                continue
    if grb_spec_index == 0:
        grb_spec_index = 2.
        logger.info('Spectrum not available. Default index set to 2.')
    else:
        logger.info('%s average spectral index: %f'%(grb_name,grb_spec_index))
    return grb_spec_index

def get_grb_position_from_site(grb_name):
    """Extract RA and DEC from Swift web site for a given GRB.

       Arguments
       ---------
       grb_name: str
               GRB name, as given (including spaces) in
               http://www.swift.ac.uk/xrt_curves/allcurvesflux.php
    """
    trigger = get_trigger_number(grb_name)
    unenh_position_url = os.path.join(PATH_UNENH_POS,trigger)
    position_page = urllib2.urlopen(unenh_position_url)
    grb_ra = 0.
    grb_dec = 0.
    ok = 0
    for line in position_page:
        if 'RA (J2000)' in line:
            if RA_DEC.search(line) is None:
                logger.info('Sorry, no position determinated for %s'\
                            %grb_name)
                ok = 1
            else:
                pass
    if ok == 0:
        position_page = urllib2.urlopen(unenh_position_url)
        for line in position_page:
            if 'RA (J2000)' in line:
                m = RA_DEC.search(line)
                grb_ra = float(m.group().replace('(',''))
            if 'Dec (J2000)' in line:
                m = RA_DEC.search(line)
                grb_dec = float(m.group().replace('(',''))
        logger.info('%s position: RA %f deg, Dec %f deg'\
                            %(grb_name,grb_ra,grb_dec))
        return grb_ra, grb_dec 
    else:
        logger.info('Set default RA = 0. and Dec = 0.')
        return grb_ra, grb_dec

def parse_light_curve_data(file_path):
    """Parse the ASCII file with the GRB light curve data from Swift.
    """

def main():
    """test the download of GRB 130427A
    """
    grb = 'GRB 130427A'
    download_swift_grb_lc_file(grb)
    download_swift_grb_lc_file(grb, light_curve='rate')


if __name__=='__main__':
    main()
