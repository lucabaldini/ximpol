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
import random

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

def get_eff_mu(_energy,mu_file_path):
    modf_en = xModulationFactor(mu_file_path)
    _modf = []
    for en in _energy:
        modf = modf_en.weighted_average([en])
        _modf.append(modf)
    _modf = numpy.array(_modf)
    return _modf

def mdp99(eff_mu, num_sig, num_bkg=0):
    """Return the MDP at the 99% confidence level.
    """
    if num_sig == 0:
        return numpy.nan
    return 4.292/(eff_mu*numpy.sqrt(num_sig + num_bkg))

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

def get_integral_flux(grb_name, delta_t=600):
    file_path = download_swift_grb_lc_file(grb_name)
    if file_path is not None:
        index = get_grb_spec_index(file_path)
        E_light_curve = parse_light_curve(file_path,num_min_data=10.)
        if E_light_curve is not None:
            scale_factor = (2. - index)/(numpy.power(MAX_ENERGY, 2. - index) - \
                                     numpy.power(MIN_ENERGY, 2. - index))
            scale_factor *= 6.242e8
            light_curve = E_light_curve.scale(scale_factor)
            t_min = light_curve.xmin()
            t_max = t_min + delta_t
            flux = light_curve.integral(t_min,t_max)
            return flux, t_min
        else:
            return None, None
    else:
        return None, None

def get_grb_mdp(grb_name, repointing=0., obs_time=100000, \
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
        E_light_curve = parse_light_curve(file_path,num_min_data=10.)
        #E_light_curve.plot(logx=True,logy=True)
        if E_light_curve is not None:
            scale_factor = (2. - index)/(numpy.power(MAX_ENERGY, 2. - index) - \
                                     numpy.power(MIN_ENERGY, 2. - index))
            scale_factor *= 6.242e8
            light_curve = E_light_curve.scale(scale_factor)
            #light_curve.plot(logx=True,logy=True)
            t_min = repointing
            if t_min < light_curve.xmin():
                logger.info('Repointing time < to the minimum time of the burst...')
                return None
            else:
                t_max = t_min + obs_time
                if t_max > light_curve.xmax():
                    t_max = light_curve.xmax()
                norm = light_curve.integral(t_min,t_max)
                _energy = numpy.linspace(MIN_ENERGY,MAX_ENERGY,1000)
                _modf = get_eff_mu(_energy,mu_file_path)
                _spectrum = get_spectrum(_energy,norm,index,aeff_file_path)
                #fmt = dict(yname='Energy Spectrum', \
                    #xname='Energy', xunits='keV', \
                    #yunits='keV$^{-1}$')
                energy_spectrum = xInterpolatedUnivariateSplineLinear(_energy,\
                                                                      _spectrum)
                #energy_spectrum.plot(logy=False, logx=False)
                num_evt = energy_spectrum.integral(MIN_ENERGY,MAX_ENERGY)
                if not math.isnan(num_evt) and num_evt != 0.:
                    logger.info('Total estimated number of events: %i'%int(num_evt))
                    mu = numpy.average(_modf,weights=_spectrum)
                    mdp = mdp99(mu,int(num_evt))
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
        
def plot_grb_mdp_vs_repoint(grb_name, _t_repoint, t_obs=100000, \
                            color='black', show=True):
    """Plot all the MDP (changing the repointing elapsed time defined in *arg) 
       for a given GRB.
    """
    mdp_list = []
    for repoint in _t_repoint:
        mdp = get_grb_mdp(grb_name,repointing=repoint,obs_time=t_obs)
        if mdp is not None:
            mdp_list.append(mdp)
        else:
            mdp_list.append(0.)
    _mdp = numpy.array(mdp_list)*100
    n_groups = len(_t_repoint)
    plt.plot(_t_repoint, _mdp, marker='o', linestyle=':', color=color)
    plt.xlabel('$t_{repoint}$ [s]')
    plt.ylabel('MDP (\%)')
    plt.title('%s -- MDP vs $t_{repoint}$, $\Delta t_{obs} =$ %i s'\
              %(grb_name,t_obs))
    if show:
        plt.show()

def plot_grb_mdp_vs_obstime(grb_name, _t_obs, t_repoint=21600, \
                            color='black', show=True):
    """Plot all the MDP (changing the repointing elapsed time defined in *arg) 
       for a given GRB.
    """
    mdp_list = []
    for obs in _t_obs:
        mdp = get_grb_mdp(grb_name,repointing=t_repoint,obs_time=obs)
        if mdp is not None:
            mdp_list.append(mdp)
        else:
            mdp_list.append(0.)
    _mdp = numpy.array(mdp_list)*100
    n_groups = len(_t_obs)
    plt.plot(_t_obs,_mdp, marker='o', linestyle=':', color=color)
    plt.xlabel('$\Delta t_{obs}$ [s]')
    plt.ylabel('MDP (\%)')
    plt.set_yscale('log')
    plt.title('%s, MDP vs $\Delta t_{obs}$, $ t_{repoint} =$ %i s'\
              %(grb_name,t_repoint))
    if show:
        plt.show()

def main():
    """Produce some plots
    """
    # If mdp_repoint = True Produce the plot of the MDP for a given GRB 
    # for different repointing times
    mdp_repoint = True
    if mdp_repoint == True:
        grb_list = ['GRB 130427A','GRB 080319B']
        plt.figure(figsize=(10, 6), dpi=80)
        for grb in grb_list:
            c = [random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]
            repointing_time = numpy.linspace(0,60000,30)
            plot_grb_mdp_vs_repoint(grb,repointing_time,show=False,color=c)
        overlay_tag()
        save_current_figure('grb_MDP',clear=False)
        #for grb in grb_list:
            #plt.figure(figsize=(10, 6), dpi=80)
            #obs_time = numpy.logspace(100,100000,20)
            #plot_grb_mdp_vs_obstime(grb,obs_time,show=False)
        plt.show()

    # If mdp_repoint = True, produce the plot of the MDP for all the Swift 
    # GRBs and a given repointing time
    all_mdp_repoint = True
    if all_mdp_repoint == True:
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
        save_current_figure('grb_MDP_repoint', clear=False)
        plt.show()
        
    # If mdp_promptflux = True, produce the plot of the correlation 
    # between MDP for all the Swift GRBs and a given repointing time
    # and the integral prompt (first 10 min) flux
    all_mdp_promptflux = True
    if all_mdp_promptflux == True:
        plt.figure(figsize=(10, 6), dpi=80)
        grb_list = get_all_swift_grb_names()        
        t_rep = 21600
        t_obs = 100000
        promt_time = 600
        mdp_list = []
        flux_list = []
        t0_list = []
        c= []
        for grb in grb_list:
            mdp = get_grb_mdp(grb,repointing=t_rep,obs_time=t_obs)
            flux, t0 = get_integral_flux(grb,delta_t=promt_time)
            if mdp is not None and flux is not None:
                mdp_list.append(mdp*100)
                flux_list.append(flux)
                t0_list.append(t0)
                if grb is 'GRB 130427A' or grb is 'GRB 080319B':
                    c.append('red')
                elif t0 > 350:
                    c.append('grey')
                else:
                    c.append('blue')
        _mdp = numpy.array(mdp_list)
        _flux = numpy.array(flux_list)
        _t0 = numpy.array(t0_list)
        plt.scatter(_mdp, _flux, s=30, marker='.', color=c)
        plt.xlabel('MDP (\%)')
        plt.ylabel('Prompt [first %i s] integral flux [keV$^{-1}$cm$^{-2}$]'\
                   %promt_time)
        plt.title('$\Delta t_{obs}=%i s,$ $t_{repoint}=%i s$'%(t_obs,t_rep))
        overlay_tag()
        save_current_figure('grb_MDP_prompt',clear=False)
        plt.show()

if __name__=='__main__':
    main()
