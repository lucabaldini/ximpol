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





IRF_NAME = 'xipe_baseline'
MIN_ENERGY = 2.
MAX_ENERGY = 10.
ENERGY_BINNING = numpy.array([MIN_ENERGY, MAX_ENERGY])

from ximpol.irf import load_arf, load_mrf
from ximpol.irf.mrf import mdp99
from ximpol.srcmodel.spectrum import int_eflux2pl_norm, xCountSpectrum

aeff = load_arf(IRF_NAME)
modf = load_mrf(IRF_NAME)




def build_count_spectrum(grb_name):
    """
    """
    file_path = download_swift_grb_lc_file(grb_name)
    if file_path is None:
        return
    logger.info('Parsing light curve data from %s...' % file_path)
    light_curve = parse_light_curve(file_path, num_min_data=10.)
    if light_curve is None:
        return
    logger.info('Done, light curve has %d point(s) between %.3f and %.3f s.' %\
                (len(light_curve.x), light_curve.xmin(), light_curve.xmax()))
    print light_curve.x, light_curve.y
    light_curve.plot(overlay=True)
    pl_index = get_grb_spec_index(file_path)
    scale_factor = int_eflux2pl_norm(1., 0.3, 10., pl_index, erg=True)
    pl_norm = light_curve.scale(scale_factor, yname='Power-law normalization',
                                yunits='keV cm$^{-2}$ s$^{-1}$')



    def energy_spectrum(E, t):
        return pl_norm(t)*numpy.power(E, -pl_index)

    return xCountSpectrum(energy_spectrum, aeff, light_curve.x)


def process_grb(grb_name, tstart=21600., duration=30000., prompt_duration=600):
    """
    """
    file_path = download_swift_grb_lc_file(grb_name)
    if file_path is None:
        return
    pl_index = get_grb_spec_index(file_path)
    light_curve = parse_light_curve(file_path, num_min_data=10.)
    if light_curve is None:
        return
    t = light_curve.x
    prompt_tstart = t[0]
    prompt_tstop = t[0] + prompt_duration
    prompt_flux = light_curve.integral(prompt_tstart, prompt_tstop)
    logger.info('Integral energy flux in %.3f--%.3f s: %.3e erg cm^{-2}' %\
                (prompt_tstart, prompt_tstop, prompt_flux))
    tstart = max(tstart, t[0])
    tstop = min(tstart + duration, t[-1])
    logger.info('Effective time interval for the MDP: %.3f--%.3f s' %\
                (tstart, tstop))
    t = t[(t >= tstart)*(t <= tstop)]
    if len(t) < 2:
        return
    
    scale_factor = int_eflux2pl_norm(1., 0.3, 10., pl_index, erg=True)
    pl_norm = light_curve.scale(scale_factor)# Fix the label.

    def energy_spectrum(E, t):
        return pl_norm(t)*numpy.power(E, -pl_index)
    
    count_spectrum = xCountSpectrum(energy_spectrum, aeff, t)
    mdp_table = count_spectrum.build_mdp_table(ENERGY_BINNING, modf)
    logger.info(mdp_table)


def process_grb_list(tstart=21600., duration=30000.):
    """
    """
    for grb_name in get_all_swift_grb_names()[:10]:
        logger.info('Processing %s...' % grb_name)
        #process_grb(grb_name, tstart, duration)
        count_spectrum = build_count_spectrum(grb_name)
        #if count_spectrum is not None:
        #    count_spectrum.plot()
    
    

def get_spectrum(_energy, norm, index):
    """Power law assumed for the energy. 
       Returns the array with the spectrum values in [KeV-1 cm-2 s-1]
       given an energy array.
    """
    return aeff(_energy)*norm*numpy.power(_energy, -index)


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

def get_grb_mdp(grb_name, repointing=21600., obs_time=100000):
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
                _energy = numpy.linspace(2.,MAX_ENERGY,1000)
                _modf = modf(_energy)
                _spectrum = get_spectrum(_energy,norm,index)
                energy_spectrum = xInterpolatedUnivariateSplineLinear(_energy,\
                                                                      _spectrum)
                #plt.xlim(1.,10.)
                #energy_spectrum.plot(logx=True,logy=True)
                
                num_evt = energy_spectrum.integral(MIN_ENERGY,MAX_ENERGY)
                if not math.isnan(num_evt) and num_evt != 0.:
                    logger.info('Total estimated number of events: %i'\
                                %int(num_evt))
                    mu = numpy.average(_modf,weights=_spectrum)
                    mdp = mdp99(mu,int(num_evt))
                    print mdp
                    if not math.isnan(mdp):
                        logger.info('%.2f--%.2f keV: %i counts (%.1e s), mu %.2f, MDP %.2f%%'\
                                    %(MIN_ENERGY,MAX_ENERGY,int(num_evt),\
                                      obs_time,mu,mdp*100))
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
                
#get_grb_mdp('GRB 130427A')

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
    plt.plot(_t_repoint, _mdp, marker='.',linestyle='-', lw=0.5, color=color,\
             label=grb_name)
    plt.xlabel('$t_{repoint}$ [s]')
    plt.ylabel('2.-10. keV MDP (%)')
    plt.title('MDP vs $t_{repoint}$, $\Delta t_{obs} =$ %i s'\
              %(t_obs))
    if show:
        plt.show()
    return _mdp, _t_repoint

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
    plt.plot(_t_obs,_mdp, marker='.',linestyle='-', lw=0.5, color=color,\
             label=grb_name)
    plt.xlabel('$\Delta t_{obs}$ [s]')
    plt.ylabel('2.-10. keV MDP (%)')
    plt.title('MDP vs $\Delta t_{obs}$, $ t_{repoint} =$ %i s'\
              %(t_repoint))
    if show:
        plt.show()
    return _mdp, _t_obs

def main():
    """Produce some plots
    """
    # If all_mdp = True, produces: 
    # 1) the plot of the MDP for all the Swift GRBs 
    #    and a given repointing time
    # 2) the plot of the correlation between MDP for all the Swift 
    #    GRBs and a given repointing time and the integral prompt 
    #    (first 10 min) flux
    all_mdp = True
    if all_mdp == True:
        grb_list = get_all_swift_grb_names()[:10]
        t_rep = 21600
        t_obs = 100000
        promt_time = 600
        mdp_list1,mdp_list2, flux_list, t0_list = [], [], [], []
        c, good_grb = [], []
        for grb in grb_list:
            mdp = get_grb_mdp(grb,repointing=t_rep,obs_time=t_obs)
            calculate_mdp(grb, t_rep, t_obs)
            flux, t0 = get_integral_flux(grb,delta_t=promt_time)
            if mdp is not None and flux is not None:
                mdp_list1.append(mdp*100)
                if t0 < 350:
                    if mdp*100 <= 15:
                        c.append('red')
                    else:
                        c.append('blue')
                    mdp_list2.append(mdp*100)
                    flux_list.append(flux)
                    t0_list.append(t0)
                else:
                    continue
        _mdp1 = numpy.array(mdp_list1)
        _mdp2 = numpy.array(mdp_list2)
        _flux = numpy.array(flux_list)
        _t0 = numpy.array(t0_list)
        # 1)------------------------------------------------------
        histo = plt.figure(figsize=(10, 6), dpi=80)
        bins = numpy.linspace(0, 100, 100)
        plt.title('%i GRBs, $\Delta t_{obs}=%i s,$ $t_{repoint}=%i s$'\
                  %(len(_mdp1),t_obs,t_rep))
        plt.hist(_mdp1, bins, alpha=0.5)
        plt.xlabel('2.-10. keV MDP (%)')
        plt.ylabel('Number of GRBs')
        overlay_tag()
        save_current_figure('all_grbs_MDP_histo', clear=False)
        plt.show()
        # 1.1)----------------------------------------------------
        histo = plt.figure(figsize=(10, 6), dpi=80)
        bins = numpy.linspace(0, 100, 100)
        plt.title('%i GRBs, $\Delta t_{obs}=%i s,$ $t_{repoint}=%i s$'\
                  %(len(_mdp1),t_obs,t_rep))
        (n, bins, patches) = plt.hist(_mdp1, bins, histtype='step', cumulative=True)
        plt.xlabel('2.-10. keV MDP (%)')
        plt.ylabel('Cumulative number of GRBs')
        for i in range(0,len(bins)):
            print 'MDP %.2f%%: %i GRBs'%(i,n[i])
        overlay_tag()
        save_current_figure('all_grbs_MDP_cumulative', clear=False)
        plt.show()
        # 2)------------------------------------------------------
        plt.figure(figsize=(10, 6), dpi=80)
        ax = plt.gca()
        plt.scatter(_mdp2, _flux, s=30, marker='.', color=c)
        plt.xlabel('2.-10. keV MDP (%)')
        plt.ylabel('[keV$^{-1}$ cm$^{-2}$]')
        plt.title('$\Delta t_{obs}=%i s,$ $t_{repoint}=%i s$'%(t_obs,t_rep))
        plt.xlim(1, 100)
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
    mdp_vs_time = False
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
            repointing_time = numpy.logspace(2,4.8,30)
            plot_grb_mdp_vs_repoint(grb,repointing_time,show=False,color=c)
        ax.legend(loc='upper left', shadow=False, fontsize='small')
        plt.plot([21600, 21600], [0, 30], 'k--', lw=1, color='green')
        plt.plot([43200, 43200], [0, 30], 'k--', lw=1,color='green')
        ax.set_yscale('log')
        ax.set_xscale('log')
        overlay_tag()
        save_current_figure('grb_MDP_vs_repoint',clear=False)
        plt.show()
        #2)------------------------------------------------------
        plt.figure(figsize=(10, 6), dpi=80)
        ax = plt.gca()
        for i,grb in enumerate(grb_list):
            obs_time = numpy.logspace(3,5,30)
            plot_grb_mdp_vs_obstime(grb,obs_time,show=False,color=color_list[i])
        ax.legend(loc='upper right', shadow=False, fontsize='small')
        plt.plot([50000, 50000], [0, 50], 'k--', lw=1, color='green')
        ax.set_yscale('log')
        ax.set_xscale('log')
        overlay_tag(x=0.5)
        save_current_figure('grb_MDP_vs_obstime',clear=False)
        plt.show()

if __name__=='__main__':
    #main()
    process_grb_list()
