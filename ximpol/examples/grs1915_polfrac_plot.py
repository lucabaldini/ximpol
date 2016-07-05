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


import os
import pyregion
import numpy

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_EXAMPLES
from ximpol import xpColor
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedMap, xBinnedModulationCube
from ximpol.core.spline import xInterpolatedUnivariateSpline
from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt


def buildspline(spindegree):
    MIN_ENERGY = 0.1
    MAX_ENERGY = 15.
    
    POL_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'ascii', 'bh_spin/Polarization_spin%s.txt'%(spindegree))
    # Build the polarization degree as a function of the energy.
    _energy, _pol_degree, _pol_angle = numpy.loadtxt(POL_FILE_PATH, unpack=True)

    #Switch to have degrees
    _pol_degree /= 100.
    
    _mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
    _energy = _energy[_mask]
    _pol_degree = _pol_degree[_mask]

    fmt = dict(xname='Energy', yname='Polarization degree')
    pol_degree_spline = xInterpolatedUnivariateSpline(_energy, _pol_degree, k=1, **fmt)

    #Pol angle in radians
    _pol_angle = numpy.deg2rad(_pol_angle)

    #Switched to have degrees and not radians
    #_pol_angle = _pol_angle

    #_mask = (_energy >= MIN_ENERGY)*(_energy <= MAX_ENERGY)
    #_energy = _energy[_mask]
    _pol_angle = _pol_angle[_mask]
    fmt = dict(xname='Energy', yname='Polarization angle [rad]')
    pol_angle_spline = xInterpolatedUnivariateSpline(_energy, _pol_angle, k=1, **fmt)

    
    return pol_degree_spline, pol_angle_spline


def fetch_mcubepath(spindegree):
    return os.path.join(XIMPOL_DATA,'grs1915_105_spin%s_mcube.fits'%spindegree)
    

#from grs1519 import SIM_DURATION
SIM_DURATION =100000

NUM_RUNS = 8


def plot(angle=False):
    spin05_pol_degree_spline, spin05_pol_angle_spline = buildspline(0.5)
    spin05_mcube =  xBinnedModulationCube(fetch_mcubepath(0.5))

    spin09_pol_degree_spline, spin09_pol_angle_spline = buildspline(0.9)
    spin09_mcube =  xBinnedModulationCube(fetch_mcubepath(0.9))

    spin998_pol_degree_spline, spin998_pol_angle_spline = buildspline(0.998)
    spin998_mcube =  xBinnedModulationCube(fetch_mcubepath(0.998))

    spin05_mcube.fit()
    spin09_mcube.fit()
    spin998_mcube.fit()
    
    spin05_fit_results = spin05_mcube.fit_results[0]
    spin09_fit_results = spin09_mcube.fit_results[0]
    spin998_fit_results = spin998_mcube.fit_results[0]

    
    plt.figure('Polarization degree')
    
    spin05_mcube.plot_polarization_degree(show=False, color='blue')
    
    spin05_pol_degree_spline.plot(color='lightblue',label='Spin 0.5', show=False)

    spin998_mcube.plot_polarization_degree(show=False, color='red')
    spin998_pol_degree_spline.plot(color='lightsalmon',label='Spin 0.998', show=False)
    
    plt.figtext(0.2, 0.85,'XIPE %s ks'%((SIM_DURATION*NUM_RUNS)/1000.),size=18)
    plt.ylim([0.00,0.045])
    plt.xlim([1,10])
    plt.legend()
    plt.show()

    if angle:
        plt.figure('Polarization angle')
    
        spin05_mcube.plot_polarization_angle(show=False, degree=True, color='blue')
        #Converting to degrees
        spin05_y = numpy.degrees(spin05_pol_angle_spline.y)
        energy = spin05_pol_angle_spline.x
        
        plt.plot(energy, spin05_y, color='lightblue',label='Spin 0.5')


        spin09_mcube.plot_polarization_angle(show=False, degree=True, color='gray')
        #Converting to degrees
        spin09_y = numpy.degrees(spin09_pol_angle_spline.y)
        energy = spin09_pol_angle_spline.x
        
        plt.plot(energy, spin09_y, color='lightgray',label='Spin 0.9')
    

        spin998_mcube.plot_polarization_angle(show=False, degree=True, color='red')
        spin998_y = numpy.degrees(spin998_pol_angle_spline.y)
        energy = spin998_pol_angle_spline.x
        
        plt.plot(energy, spin998_y, color='lightsalmon',label='Spin 0.998')
        #spin998_pol_angle_spline.plot(color='lightsalmon',label='Spin 0.998', show=False)
    
        plt.figtext(0.2, 0.85,'XIPE %s ks'%((SIM_DURATION*NUM_RUNS)/1000.),size=18)
        plt.xlim([1,10])
        plt.ylim([40,200])
        plt.legend()
        plt.show()

def plotmdp():
    spin00_pol_degree_spline = buildspline(0.5)
    spin00_mcube =  xBinnedModulationCube(fetch_mcubepath(0.5))

    spin998_pol_degree_spline = buildspline(0.998)
    spin998_mcube =  xBinnedModulationCube(fetch_mcubepath(0.998))

    spin00_mcube.fit()
    spin998_mcube.fit()
    
    spin00_fit_results = spin00_mcube.fit_results[0]
    spin998_fit_results = spin998_mcube.fit_results[0]

    
    plt.figure('MDP')
    
    spin00_mdp = spin00_mcube.mdp99[:-1]
    spin998_mdp = spin998_mcube.mdp99[:-1]
    emean = spin00_mcube.emean[:-1]
    emin =  spin00_mcube.emin[:-1]
    emax =  spin00_mcube.emax[:-1]
    width = (emax-emin)/2.
    plt.errorbar(emean,spin00_mdp,xerr=width, label='Spin 0.5',marker='o',linestyle='--')
    
    plt.errorbar(emean,spin998_mdp,xerr=width, label='Spin 0.998',marker='o',linestyle='--')
            
    plt.figtext(0.2, 0.85,'XIPE %s ks'%((SIM_DURATION*NUM_RUNS)/1000.),size=18)
    plt.xlim([1,10])
    
    plt.ylabel('MPD 99\%')
    plt.xlabel('Energy (keV)')
    plt.legend()
    plt.show()



    
if __name__ == '__main__':
    plot(True)
    plotmdp()
