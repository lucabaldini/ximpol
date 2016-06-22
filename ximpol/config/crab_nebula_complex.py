#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
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


import numpy
import os

from ximpol.srcmodel.roi import xExtendedSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import xPolarizationMap, constant
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.logging_ import logger
from ximpol import XIMPOL_CONFIG


"""Configuration file for a semi-realistic model of Cas A.

The spectral model is taken from
E.A. Helder and J. Vink, "Characterizing the non-thermal emission of Cas A",
Astrophys.J. 686 (2008) 1094--1102, http://arxiv.org/abs/0806.3748,
which seems to be one of the few instances where an actual spectrum in physical
units (i.e., not a count spectrum) is presented.
We grabbed by hand the data points and we're calling "thermal" whatever is in
the lines and "non-thermal" the rest.

We have two images of Cas A, at low (1.5--3.0 keV) and high (4.0--6.0 keV) and,
due to the absence of lines between 4 and 6 keV we're attaching the latter to
the non-thermal spectrum and the former to the thermal component.

The polarization map is a simple geometrical, radially-symmetric, model.
"""

import numpy
import os

from ximpol.srcmodel.roi import xExtendedSource, xROIModel
from ximpol.srcmodel.spectrum import power_law
from ximpol.srcmodel.polarization import constant
from ximpol import XIMPOL_CONFIG

# Location of the ROI:
ROI_MODEL = xROIModel(83.633083, 22.014500)
# Image:
img_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'crab_0p3_10p0_keV.fits')
# Spectrum

# The idea is to have multiple extended sources, one per region.
# For example to create multiple region polarization map, from the bin diurectory:
#./xppolmap.py --image ../config/fits/crab_0p3_10p0_keV.fits --region ../config/fits/crab_nebula_complex.reg --output ../config/fits/crab_nebula_complex.fits
##############################

base_names=['crab_nebula_pmax050_reg000','crab_nebula_pmax050_reg001','crab_nebula_pmax050_reg002']
#base_names=['crab_nebula_pmax050_reg000']
polarization_maps=[]
for i,base_name in enumerate(base_names):
    energy_spectrum = power_law(9.59/len(base_names), 2.108) # This ineeds to be checked!
    
    pol_mapx_path = os.path.join(XIMPOL_CONFIG, 'fits', '%s_x.fits' % base_name)
    pol_mapy_path = os.path.join(XIMPOL_CONFIG, 'fits', '%s_y.fits' % base_name)
    polarization_map = xPolarizationMap(pol_mapx_path, pol_mapy_path)
    polarization_maps.append(polarization_map)
    def nebula_polarization_angle(E, t, ra, dec):
        return polarization_map.polarization_angle(ra, dec)

    def nebula_polarization_degree(E, t, ra, dec):
        return polarization_map.polarization_degree(ra, dec)
    
    _source=xExtendedSource('Crab Nebula (%d)' % i,
                            img_file_path,
                            energy_spectrum,
                            nebula_polarization_degree,
                            nebula_polarization_angle)
    ROI_MODEL.add_source(_source)
    pass



if __name__ == '__main__':
    print(ROI_MODEL)

def display():
    """Display the source model.
    """
    from ximpol.utils.matplotlib_ import pyplot as plt
    from ximpol.srcmodel.img import xFITSImage

    print(ROI_MODEL)
    #fig = plt.figure('Energy spectrum')
    #for src in  ROI_MODEL.values():
        #src.energy_spectrum.plot(logy=True, show=False, label=src.name) 
    #plt.legend(bbox_to_anchor=(0.95, 0.95))
    #    fig = thermal_component.image.plot(show=False)
    #xFITSImage.add_label(fig, 'Chandra 1.5-3.0 keV')
    #fig = nonthermal_component.image.plot(show=False)
    #xFITSImage.add_label(fig, 'Chandra 4.0-6.0 keV')
    img = xFITSImage(img_file_path)
    fig = img.plot(show=False)
    for polarization_map in polarization_maps:
        polarization_map.build_grid_sample(ROI_MODEL.ra, ROI_MODEL.dec,num_points=100)
        polarization_map.overlay_arrows(fig)
    plt.show()


if __name__ == '__main__':
    display()
