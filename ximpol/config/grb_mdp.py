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

import os
import sys
import numpy
import math

from ximpol import XIMPOL_IRF
from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.mrf import xModulationFactor
from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.config.grb_swift_download import download_swift_grb_lc_file
from ximpol.config.grb_swift_download import get_all_swift_grb_names
from ximpol.config.grb_utils import parse_light_curve
from ximpol.config.grb_utils import get_grb_spec_index

MAX_ENERGY = 10.
MIN_ENERGY = 0.3

def get_eff_mu(emin,emax,mu_file_path):
    modf_en = xModulationFactor(mu_file_path)
    total_modf = modf_en.integral(emin,emax)/(emax-emin)
    return total_modf

def mdp99(eff_mu, num_sig, num_bkg=5.):
    """Return the MDP at the 99% confidence level.
    """
    if num_sig == 0:
        return numpy.nan
    return 4.292/(eff_mu*numpy.sqrt(num_sig + num_bkg))#/num_sig

def get_spectrum(_energy, norm, index, aeff_file_path):
    """Power law assumed for the energy. 
       Returns the array with the spectrum values in [KeV-1 cm-2 s-1]
       given an energy array.
    """
    aeff_spline = xEffectiveArea(aeff_file_path)
    #getting the effective area as a function of E for given off-axis angle
    aeff_list = []
    _aeff = aeff_spline.eval_(_energy,theta=0.)
    return _aeff*norm*numpy.power(_energy, -index)

def get_grb_mdp(grb_name, repointing=21600, obs_time=100000, \
                aeff_file_path=os.path.join(XIMPOL_IRF,'fits',\
                                            'xipe_baseline.arf'), \
                mu_file_path=os.path.join(XIMPOL_IRF,'fits',\
                                             'xipe_baseline.mrf')):
    """Calculate the MDP, given GRB name, the repointing elapsed time [sec]
       after the trigger, and the oservationrange of time [sec]
    """
    file_path = download_swift_grb_lc_file(grb_name)
    if file_path is not None:
        index = get_grb_spec_index(file_path)
        E_light_curve = parse_light_curve(file_path,num_min_data=10)
        if E_light_curve is not None:
            scale_factor = (2. - index)/(numpy.power(MAX_ENERGY, 2. - index) - \
                                     numpy.power(MIN_ENERGY, 2. - index))
            scale_factor *= 6.242e8
            light_curve = E_light_curve.scale(scale_factor)
            t_min = repointing
            t_max = repointing + obs_time
            if t_max > light_curve.xmax():
                t_max = light_curve.xmax()
            norm = light_curve.integral(t_min,t_max)
            _energy = numpy.linspace(MIN_ENERGY,MAX_ENERGY,1000)
            _spectrum = get_spectrum(_energy,norm,index,aeff_file_path)
            energy_spectrum = xInterpolatedUnivariateSplineLinear(_energy,\
                                                                  _spectrum)
            num_evt = energy_spectrum.integral(MIN_ENERGY,MAX_ENERGY)
            if not math.isnan(num_evt):
                logger.info('Total estimated number of events: %i'%int(num_evt))
                modf = get_eff_mu(MIN_ENERGY,MAX_ENERGY,mu_file_path)
                mdp = mdp99(modf,int(num_evt))
                if not math.isnan(mdp):
                    return mdp
                else:
                    return None
            else:
                logger.info('Total number of events cannot be estimated...')
                return None
        else:
            return None
    else:
        return None
        
def plot_grb_mdp(grb_name,_t,show=True):
    """Plot all the MDP (changing the repointing elapsed time defined in *arg) 
       for a given GRB.
    """
    mpd_list = []
    for repoint in _t:
        mdp = get_grb_mdp(grb_name,repointing=repoint)
        if mdp is not None:
            mpd_list.append(mdp)
        else:
            mdp_list.append(0.)
    _mdp = numpy.array(mpd_list)
    n_groups = len(_t)
    rects1 = plt.bar(_t, _mdp*100, 10000.,
                     alpha = 0.5,
                     color='r')
    plt.xlabel('Seconds after the trigger')
    plt.ylabel('MDP (\%)')
    plt.title('MDP for %s' %grb_name)
    plt.xticks(_t + 5000.)
    overlay_tag(x=0.60)
    save_current_figure('grb_MDP_%s'%grb_name.replace(' ','').replace('.','-'),
                        clear=False)
    if show:
        plt.show()

def main():
    """Produce some plots
    """
    #Produce the plot of the MDP for a given GRB 
    #for different repointing times
    grb_name = 'GRB 130427A'
    repointing_time = numpy.linspace(21600,86400,3)
    plot_grb_mdp(grb_name,repointing_time,show=False)
    """
    #Produce the plot of the MDP for all the Swift 
    #GRBs and a given repointing time
    histo = plt.figure(figsize=(10, 6), dpi=80)
    grb_list = get_all_swift_grb_names()
    mdp_list = []
    bins = numpy.linspace(0, 100, 100)
    t = 21600
    for grb in grb_list:
        mdp = get_grb_mdp(grb,repointing=t)
        if mdp is not None:
            mdp_list.append(mdp*100)
    _mdp = numpy.array(mdp_list)
    plt.title('Sample of %i GRBs, MDPs if repointing after %i seconds'\
              %(len(_mdp),t))
    plt.hist(_mdp, bins, alpha=0.5)
    plt.xlabel('MDP (\%)')
    plt.ylabel('Number of GRBs')
    overlay_tag()
    save_current_figure('grb_MDP_%s'%grb_name.replace(' ','').replace('.','-'),
                        clear=False)"""
    plt.show()

if __name__=='__main__':
    main()
