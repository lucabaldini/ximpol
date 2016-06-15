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

def get_integral_flux(grb_name,delta_t=600):
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
            return flux, light_curve.xmin()
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
                energy_spectrum = xInterpolatedUnivariateSplineLinear(_energy,\
                                                                      _spectrum)
                num_evt = energy_spectrum.integral(MIN_ENERGY,MAX_ENERGY)
                if not math.isnan(num_evt) and num_evt != 0.:
                    logger.info('Total estimated number of events: %i'\
                                %int(num_evt))
                    mu = numpy.average(_modf,weights=_spectrum)
                    mdp = mdp99(mu,int(num_evt))
                    if not math.isnan(mdp):
                        logger.info('%lf--%lf keV: %i counts in %i s, mu %lf, MDP %lf'\
                                    %(MIN_ENERGY,MAX_ENERGY,int(num_evt),obs_time,mu,mdp))
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
    plt.plot(_t_repoint, _mdp, marker='.',linestyle='--', lw=0.5, color=color,\
             label=grb_name)
    plt.xlabel('$t_{repoint}$ [s]')
    plt.ylabel('MDP (\%)')
    plt.title('MDP vs $t_{repoint}$, $\Delta t_{obs} =$ %i s'\
              %(t_obs))
    if show:
        plt.show()

target = []
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
    target = [i for i in _mdp if i<0.09*_mdp[0]]
    plt.plot(_t_obs,_mdp, marker='.',linestyle='--', lw=0.5, color=color,\
             label=grb_name)
    plt.xlabel('$\Delta t_{obs}$ [s]')
    plt.ylabel('MDP (\%)')
    plt.title('MDP vs $\Delta t_{obs}$, $ t_{repoint} =$ %i s'\
              %(t_repoint))
    print _mdp
    if show:
        plt.show()

def main():
    """Produce some plots
    """
    # If all_mdp = True, produces: 
    # 1) the plot of the MDP for all the Swift GRBs 
    #    and a given repointing time
    # 2) the plot of the correlation between MDP for all the Swift 
    #    GRBs and a given repointing time and the integral prompt 
    #    (first 10 min) flux
    good_grb = []
    all_mdp = False
    if all_mdp == True:
        grb_list = get_all_swift_grb_names()
        bins = numpy.linspace(0, 100, 100)
        t_rep = 21600
        t_obs = 100000
        promt_time = 600
        mdp_list, flux_list, t0_list = [], [], []
        c, good_grb = [], []
        for grb in grb_list:
            mdp = get_grb_mdp(grb,repointing=t_rep,obs_time=t_obs)
            flux, t0 = get_integral_flux(grb,delta_t=promt_time)
            if mdp is not None and flux is not None:
                mdp_list.append(mdp*100)
                flux_list.append(flux)
                t0_list.append(t0)
                if mdp*100 <= 15:
                    c.append('red')
                    good_grb.append(grb)
                elif t0 > 350:
                    c.append('grey')
                else:
                    c.append('blue')
        _mdp = numpy.array(mdp_list)
        _flux = numpy.array(flux_list)
        _t0 = numpy.array(t0_list)
        # 1)------------------------------------------------------
        histo = plt.figure(figsize=(10, 6), dpi=80)
        plt.title('Sample of %i GRBs, MDPs if repointing after %i seconds'\
                  %(len(_mdp),t_rep))
        plt.hist(_mdp, bins, alpha=0.5)
        plt.xlabel('MDP (\%)')
        plt.ylabel('Number of GRBs')
        overlay_tag()
        save_current_figure('grb_MDP_repoint', clear=False)
        # 2)------------------------------------------------------
        plt.figure(figsize=(10, 6), dpi=80)
        ax = plt.gca()
        plt.scatter(_mdp, _flux, s=30, marker='.', color=c)
        plt.xlabel('MDP (\%)')
        plt.ylabel('Prompt [first %i s] integral flux [keV$^{-1}$cm$^{-2}$]'\
                   %(promt_time,21600))
        plt.title('$\Delta t_{obs}=%i s,$ $t_{repoint}=%i s$'%(t_obs,t_rep))
        ax.set_yscale('log')
        ax.set_xscale('log')
        overlay_tag()
        save_current_figure('grb_MDP_prompt',clear=False)
        plt.show()
    
    # If mdp_vs_time = True Produces: 
    # 1) the plot of the MDP for a given GRB 
    #    as a function of the repointing time
    # 2) the plot of the MDP for a given GRB 
    #    as a function of the observation duration
    mdp_vs_time = True
    color_list = []
    if mdp_vs_time == True:
        grb_list = ['GRB 060729', 'GRB 080411', 'GRB 091127', 'GRB 111209A',\
                    'GRB 120711A', 'GRB 130427A', 'GRB 130505A', 'GRB 130907A',\
                    'GRB 150403A']
        #1)------------------------------------------------------
        plt.figure(figsize=(10, 6), dpi=80)
        ax = plt.gca()
        for grb in grb_list:
            c = [random.uniform(0,1),random.uniform(0,1),random.uniform(0,1)]
            color_list.append(c)
            repointing_time = numpy.linspace(100,60000,30)
            plot_grb_mdp_vs_repoint(grb,repointing_time,show=False,color=c)
        ax.legend(loc='upper left', shadow=False)
        plt.plot([21600, 21600], [0, 30], 'k--', lw=1, color='green')
        plt.plot([43200, 43200], [0, 30], 'k--', lw=1,color='green')
        overlay_tag()
        save_current_figure('grb_MDP_vs_repoint',clear=False)
        #2)------------------------------------------------------
        plt.figure(figsize=(10, 6), dpi=80)
        ax = plt.gca()
        for i,grb in enumerate(grb_list):
            obs_time = numpy.linspace(1000,50000,20)
            plot_grb_mdp_vs_obstime(grb,obs_time,show=False,color=color_list[i])
        ax.legend(loc='upper right', shadow=False)
        overlay_tag(x=0.5)
        save_current_figure('grb_MDP_vs_obstime',clear=False)
        plt.show()

if __name__=='__main__':
    main()
