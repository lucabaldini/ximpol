#!/urs/bin/env python
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


import sys
import os
import numpy as np

from ximpol.srcmodel.gabs import get_column_density, get_galactic_absorption
from ximpol.utils.logging_ import logger,abort
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


def XsectionISM_array(energy_array):
    """Returns x-sections array corresponding to a given energy array
    """
    from ximpol.srcmodel.gabs import XsectionISM
    xsec_array = []
    for en in energy_array:
        xsec = XsectionISM(en)
        xsec_array.append(xsec)
    xsec_array = np.array(xsec_array)
    return xsec_array 


def XsectionISM_energy_plot(energy_array, show=True):
    """Plots the x-sections*E^3 in 10^{-24}cm^2 as a function of the 
       energy in keV
    """
    xsec_array = XsectionISM_array(energy_array)
    xsec_array_e3 = xsec_array*(energy_array**3)*10**24
    title = 'X-Ray Absorption Cross-section in the ISM'
    plt.figure(title)
    plt.plot(energy_array,xsec_array_e3)
    plt.title(title)
    plt.xlabel('E [KeV]')
    plt.xscale('log')
    plt.axis([0.03,10.,0.,1000])
    plt.ylabel('$\sigma \cdot$ E$^3$ [10$^{-24}$ cm$^2$]')
    if show:
        plt.show()

def TransmissionFactor_array(energy_array, column_density):
    """Returns the transmission factors array corresponding to a 
       given energy array 
    """
    trans_factor_array = []
    xsec_array = XsectionISM_array(energy_array)
    for energy in energy_array:
        trans_factor = get_galactic_absorption(energy,column_density)
        trans_factor_array.append(trans_factor)
    trans_factor_array = np.array(trans_factor_array)
    return trans_factor_array


def TransmissionFactor_energy_plot(energy_array, column_density, show=True):
    """Plots the transmission factor as a function of the 
       energy in keV  
    """
    trans_factor_array = TransmissionFactor_array(energy_array,column_density)
    title = 'X-Ray Transmission Factor for Nh=%e'%column_density
    plt.figure(title)
    plt.plot(energy_array,trans_factor_array)
    plt.title(title)
    plt.xlabel('E [KeV]')
    plt.xscale('log')
    plt.axis([0.03,10.,0.,1.])
    plt.ylabel('Transmission Factor')
    if show:
        plt.show()


def main(interactive=False):
    """
    """
    RA,DEC = 10.684, 41.269
    column_density = get_column_density(RA,DEC,'LAB')
    E_min, E_max = 0.03, 10.
    logger.info('Energy values [keV] range between %f and %f'%(E_min,E_max))
    E_bins = np.linspace(E_min,E_max,1000)
    logger.info('Plotting the X-ray absorption cross-section in the ISM...')
    XsectionISM_energy_plot(E_bins,show=False)
    overlay_tag()
    save_current_figure('absorption_cross_section.png', clear=False)
    logger.info('Plotting the X-ray transmission factor as a function of the energy...')
    TransmissionFactor_energy_plot(E_bins,column_density,show=False)
    overlay_tag(x=0.55)
    save_current_figure('test_transmission_factor_NH%d.png'%column_density,
                        clear=False)
    if interactive:
        plt.show()

        
if __name__ == '__main__':
    main(interactive=sys.flags.interactive)
