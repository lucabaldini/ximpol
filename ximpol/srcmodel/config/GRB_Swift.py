#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
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
import sys,os
from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.srcmodel.spectrum import power_law, constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_SRCMODEL

roi_ra,roi_dec=173.13, 27.71
grb_ra,grb_dec=173.1362, 27.7129

ROI_MODEL = xROIModel(roi_ra,roi_dec)
import scipy
import scipy.signal
def parse(ascii_file_name,scale=3):
    time_array = []
    flux_array = []
    for l in file(ascii_file_name,'r'):
        try:
            t,tp,tm,f,fp,fm=[float(item) for item in l.split()]
            time_array.append(t)
            flux_array.append(f)
        except:
            pass
        pass
    time_array=numpy.array(time_array)
    flux_array=numpy.array(flux_array)

    flux_array=scipy.signal.medfilt(flux_array,scale)
    return time_array,flux_array
    



time_array,flux_array=parse(os.path.join(XIMPOL_SRCMODEL,'ascii/GRB130427_Swift.dat'),scale=11)

source1 = xPointSource(name='GRB', ra=grb_ra, dec=grb_dec)
integral_flux    = xInterpolatedUnivariateSplineLinear(time_array,flux_array,'Time','s','Flux','erg/cm^2/s')
erg2kev=6.242e+8
def gamma(t):
    return 1.66

def C(t,emin=0.3,emax=10.):
    return (gamma(t)-1.)*integral_flux(t)/(numpy.power(emin,-gamma(t)+1.0)-numpy.power(emax,-gamma(t)+1.0))*erg2kev

def dNde(e,t):
    return C(t)*numpy.power(e,-gamma(t))

def polarization_degree(e,t):
    return 0.6*(1.0-t/time_array[-1])

source1.spectrum = dNde
source1.polarization_degree = polarization_degree
source1.polarization_angle  = constant(0.0)
ROI_MODEL.add_source(source1)
if __name__=='__main__':
    from matplotlib import pyplot as plt
    fig=plt.figure(figsize=(10,10))
    fig.add_subplot(2,1,1)
    plt.plot(time_array,dNde(1.,time_array))
    plt.plot(time_array,dNde(2.,time_array))
    plt.xscale('log')
    plt.yscale('log')
    fig.add_subplot(2,1,2)
    plt.plot(time_array,polarization_degree(1.,time_array))
    plt.xscale('log')
    plt.show()
