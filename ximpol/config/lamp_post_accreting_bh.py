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


import numpy
import os

from ximpol.srcmodel.roi import xPointSource, xROIModel
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol import XIMPOL_CONFIG


def parse(ascii_file_name, eMin = 1, eMax = 80, scale=3):
    ene_array = []
    flux_array = []
    pol_degree_array = []
    pol_angle_array = []
    for l in file(ascii_file_name,'r'):
        try:
            E,flux,pol_deg,pol_angle=[float(item) for item in l.split()]
            if E >= eMin and E <= eMax:
                ene_array.append(E)
                flux_array.append(flux)
                pol_degree_array.append(pol_deg)
                pol_angle_array.append(pol_angle)
        except:
            pass
        pass
    ene_array = numpy.array(ene_array)
    flux_array = numpy.array(flux_array)
    pol_degree_array = numpy.array(pol_degree_array)
    pol_angle_array = numpy.array(pol_angle_array)
    return ene_array, flux_array, pol_degree_array, pol_angle_array


"""
In the so called 'lamp-post' model, the unpolarized emission from a primary
source illuminates the accretion disk of a black hole (BH), where it is
reprocessed. Part of it is emitted toward the observer, with a polarization
degree and angle depending on the exact geometry of the model.

Main parameters of the 'lamp post' emission model are:
    - the height of the primary source above to the accretion disk, in units
      of the gravitational radius
    - the observer inclination angle
    - the angular momentum per unit mass of the black hole in units of
      the gravitational radius. Spin == 0 correspond to a static
      BH (Schwarzchild metric), Spin == 1 to an extremely rotating one
      (Kerr metric).

Here we are using spectrum, polarization degree and polarization angle
for a model of the Seyfert 1 galaxy MCG-6-30-15 presented in
Dovciak et al. 'Light bending scenario for accreting black holes in x-ray
polarimetry' (2011), provided by Fabio Muleri.

The required info are stored as columns in several .dat files.
Each file corresponds to a different choice of the main input parameters.
"""

# spin = 0, h =  3, inclin_angle = 30
data_file_name = 'lamp_pol_0_30_003_reduced.dat'
data_file_path = os.path.join (XIMPOL_CONFIG, 'ascii', data_file_name)

eMin = 1.
eMax = 100.
ene_array, flux_array, pol_degree_array, pol_angle_array = \
                            parse(data_file_path, eMin, eMax)

flux_de = xInterpolatedUnivariateSplineLinear(ene_array, flux_array, \
                                              'Energy', 'keV',\
                                              'Flux', 'ph/cm^2/s')
pol_degree_de = xInterpolatedUnivariateSplineLinear(ene_array,\
                                                    pol_degree_array,
                                                    'Energy', 'keV',\
                                                    'Polarization degree')
pol_angle_de = xInterpolatedUnivariateSplineLinear(ene_array, pol_angle_array,\
                                                   'Energy', 'keV',\
                                                   'Polarization angle',\
                                                   'degree')

#convertion factor from arbitrary units to ph/cm2/s (provided by Fabio Muleri)
flux_normalization_factor = 6.0067e-2
# NOTE: Not sure these are the correct input units, however, since nothing
# in this model is time-dependent, we can temporarily change the observation
# time to get every desired normalization.

def dNde (e,t):
    return  flux_normalization_factor*flux_de(e)

def polarization_degree_de (e,t, ra, dec):
    return pol_degree_de (e)

def polarization_angle_de (e,t, ra, dec):
    return numpy.radians(pol_angle_de(e)) #the code works with radians

ra_MCG  = 203.973
dec_MCG = -34.2960
ROI_MODEL = xROIModel(ra_MCG, dec_MCG)

lamp_post_source = xPointSource('MCG-6-30-15', ra_MCG, dec_MCG, dNde,\
                                 polarization_degree_de, polarization_angle_de)

ROI_MODEL.add_source(lamp_post_source)

if __name__=='__main__':
    from matplotlib import pyplot as plt
    x_array = ene_array[0:550]
    fig=plt.figure(figsize=(10,10))
    fig.add_subplot(3,1,1)
    plt.plot( x_array, dNde( x_array, 1. ) )
    plt.xscale('log')
    plt.yscale('log')
    plt.ylabel('Flux (ph/cm^2/s)')
    fig.add_subplot(3,1,2)
    plt.plot( x_array, polarization_degree_de( x_array, 1., 1., 1. ) )
    plt.xscale('log')
    plt.ylabel('Polarization degree')
    fig.add_subplot(3,1,3)
    plt.plot( x_array, \
                  numpy.degrees( polarization_angle_de( x_array, 1.,1.,1. ) ) )
    plt.xscale('log')
    plt.xlabel('Energy (keV)')
    plt.ylabel('Polarization angle')
    plt.tight_layout()
    plt.show()
