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
from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt


from ximpol.config.cyg_x1_wedge_corona import pol_degree_spline as pol_degree_spline_wedge

from ximpol.config.cyg_x1_wedge_corona import pol_angle_spline as  pol_angle_spline_wedge


from ximpol.config.cyg_x1_spherical_corona import pol_degree_spline as pol_degree_spline_spherical

from ximpol.config.cyg_x1_spherical_corona import pol_angle_spline as  pol_angle_spline_spherical

from cyg_x1 import SIM_DURATION

NUM_RUNS = 3
#SPHERICAL_CFG_FILE = os.path.join(XIMPOL_CONFIG, 'cyg_x1_spherical.py')


WEDGE_MCUBE_PATH = os.path.join(XIMPOL_DATA,'cyg_x1_wedge_corona_40_inclination_mcube.fits')

SPHERICAL_MCUBE_PATH = os.path.join(XIMPOL_DATA,'cyg_x1_spherical_corona_40_inclination_mcube.fits')


def plot():
    
    spherical_mcube =  xBinnedModulationCube(SPHERICAL_MCUBE_PATH)
    wedge_mcube =  xBinnedModulationCube(WEDGE_MCUBE_PATH)

    spherical_mcube.fit()
    wedge_mcube.fit()
    
    spherical_fit_results = spherical_mcube.fit_results[0]
    wedge_fit_results = wedge_mcube.fit_results[0]

    
    plt.figure('Polarization degree')
    
    spherical_mcube.plot_polarization_degree(show=False, color='blue')
    pol_degree_spline_spherical.plot(color='lightblue',label='Spherical corona model (40 degree inclination)', show=False)

    wedge_mcube.plot_polarization_degree(show=False, color='red')
    pol_degree_spline_wedge.plot(color='lightsalmon',label='Wedge corona model  (40 degree inclination)', show=False)
    
    plt.figtext(0.2, 0.85,'XIPE %s ks'%((SIM_DURATION*NUM_RUNS)/1000.),size=18)
    plt.xlim([1,10])
    plt.legend()
    plt.show()


    
if __name__ == '__main__':
    plot()
