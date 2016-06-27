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
from astropy.io import fits

from ximpol import XIMPOL_IRF
from ximpol import XIMPOL_DATA
from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.mrf import xModulationFactor
from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.config.grb_swift_download import download_swift_grb_lc_file
from ximpol.config.grb_swift_download import get_all_swift_grb_names
from ximpol.config.grb_utils import parse_light_curve
from ximpol.config.grb_utils import get_grb_spec_index, get_grb_position
from ximpol.core.fitsio import xBinTableHDUBase, xPrimaryHDU


IRF_NAME = 'xipe_baseline'
MIN_ENERGY = 2.
MAX_ENERGY = 10.
ENERGY_BINNING = numpy.array([MIN_ENERGY, MAX_ENERGY])
OUTFILE = os.path.join(XIMPOL_DATA,'GRBmainInfos.fits')

from ximpol.irf import load_arf, load_mrf
from ximpol.irf.mrf import mdp99
from ximpol.srcmodel.spectrum import int_eflux2pl_norm, xCountSpectrum

aeff = load_arf(IRF_NAME)
modf = load_mrf(IRF_NAME)


class xBinTableGRBmain(xBinTableHDUBase):

    NAME = 'GRB_MAIN'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('NAME'        , 'A20', None         , 'grb name'),
        ('ENERGY_LO'   , 'E',    'keV'       , 'energy low'),
        ('ENERGY_HI'   , 'E',    'keV'       , 'energy high'),
        ('RA'          , 'E',   'degrees'    , 'grb right ascension'),
        ('DEC'         , 'E',   'degrees'    , 'grb declination'),
        ('INDEX'       , 'E',   None         , 'late spectral index'),
        ('START'       , 'D',   's'          , 'observation start time'),
        ('STOP'        , 'D',   's'          , 'observation stop time'),
        ('PROMPT_FLUX' , 'E',   'erg/cm2/s'  , 'grb prompt flux'),
        ('PROMPT_START', 'D',   's'          , 'grb prompt flux'),
        ('PROMPT_STOP' , 'D',   's'          , 'grb prompt flux'),
        ('GRB_START'   , 'D',   's'          , 'grb start time'),
        ('GRB_STOP'    , 'D',   's'          , 'grb stop time'),
        ('EFFECTIVE_MU', 'E',    None        , 'effective modulation factor'),
        ('COUNTS'      , 'J',    None        , 'total counts'),
        ('MDP 99%'     , 'E',    None        , 'mdp')
    ]

def build_grb_fits_file(data,outfile):
    primary_hdu = xPrimaryHDU()
    grb_info_hdu = xBinTableGRBmain(data)
    hdu_list = fits.HDUList([primary_hdu, grb_info_hdu])
    hdu_list.info()
    logger.info('Writing GRB main infos table to %s...' % outfile)
    hdu_list.writeto(outfile, clobber=True)
    logger.info('Done.')

def process_grb(grb_name, tstart=21600., duration=30000., prompt_duration=600):
    """
    """
    file_path = download_swift_grb_lc_file(grb_name)
    if file_path is None:
        return None
    pl_index = get_grb_spec_index(file_path)
    ra, dec = get_grb_position(file_path)
    light_curve = parse_light_curve(file_path, num_min_data=5.)
    if light_curve is None:
        return None
    t = light_curve.x
    grb_start, prompt_tstart = t[0], t[0]
    grb_stop = t[-1]
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
        return None
    scale_factor = int_eflux2pl_norm(1., 0.3, 10., pl_index, erg=True)
    pl_norm = light_curve.scale(scale_factor)# Fix the label.
    def energy_spectrum(E, t):
        return pl_norm(t)*numpy.power(E, -pl_index)
    count_spectrum = xCountSpectrum(energy_spectrum, aeff, t)
    mdp_table = count_spectrum.build_mdp_table(ENERGY_BINNING, modf)
    logger.info(mdp_table)
    mdp = mdp_table.mdp_values()[0]
    eff_mu = [row.mu_effective for row in mdp_table.rows]
    counts = [row.num_signal for row in mdp_table.rows]
    grb_values = [ra, dec, pl_index, tstart, tstop, prompt_flux, prompt_tstart,\
                  prompt_tstop, grb_start, grb_stop, eff_mu[0], counts[0], mdp]
    return grb_values

def process_grb_list(tstart=21600., duration=30000., prompt_duration=600):
    """
    """
    name = numpy.array([], dtype=str)
    e_low = numpy.array([])
    e_high = numpy.array([])
    ra = numpy.array([])
    dec = numpy.array([])
    index = numpy.array([])
    start = numpy.array([])
    stop = numpy.array([])
    prompt_start = numpy.array([])
    prompt_stop = numpy.array([])
    prompt_flux = numpy.array([])
    grb_start = numpy.array([])
    grb_stop = numpy.array([])
    eff_mu = numpy.array([])
    counts = numpy.array([], dtype=numpy.int64)
    mdp = numpy.array([])
    for grb_name in get_all_swift_grb_names():
        logger.info('Processing %s...' % grb_name)
        grb_values = process_grb(grb_name,tstart=tstart,duration=duration,\
                                 prompt_duration=prompt_duration)
        if grb_values is None:
            continue
        name = numpy.append(name,[grb_name])
        e_low = numpy.append(e_low,[MIN_ENERGY])
        e_high = numpy.append(e_high,[MAX_ENERGY])
        ra = numpy.append(ra,[grb_values[0]])
        dec = numpy.append(dec,[grb_values[1]])
        index = numpy.append(index,[grb_values[2]])
        start = numpy.append(start,[grb_values[3]])
        stop = numpy.append(stop,[grb_values[4]])
        prompt_flux = numpy.append(prompt_flux,[grb_values[5]])
        prompt_start = numpy.append(prompt_start,[grb_values[6]])
        prompt_stop = numpy.append(prompt_stop,[grb_values[7]])
        grb_start = numpy.append(grb_start,[grb_values[8]])
        grb_stop = numpy.append(grb_stop,[grb_values[9]])
        eff_mu = numpy.append(eff_mu,[grb_values[10]])
        counts = numpy.append(counts,[grb_values[11]])
        mdp = numpy.append(mdp,[grb_values[12]])
    grb_values_array = [name,e_low,e_high,ra,dec,index,start,stop,prompt_flux,\
                        prompt_start,prompt_stop,grb_start,grb_stop,eff_mu,\
                        counts,mdp]
    return grb_values_array
        

def get_spectrum(_energy, norm, index):
    """Power law assumed for the energy.
       Returns the array with the spectrum values in [KeV-1 cm-2 s-1]
       given an energy array.
    """
    return aeff(_energy)*norm*numpy.power(_energy, -index)


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
    # If process_grb_mdp = True, produces a fits file with all the 
    # main infos on each grb
    process_grb_mdp = True
    if process_grb_mdp == True:
        data = process_grb_list()
        build_grb_fits_file(data,OUTFILE)

    # 1) the plot of the MDP for all the Swift GRBs
    #    and a given repointing time
    # 2) the cumulative of the previous histogram
    # 3) the plot of the correlation between MDP for all the Swift
    #    GRBs and a given repointing time and the integral prompt
    #    (first 10 min) flux

    # 1)------------------------------------------------------
    plt.figure(figsize=(10, 6), dpi=80)
    bins = numpy.linspace(0, 100, 100)
    hdulist = fits.open(OUTFILE)
    grbdata = hdulist[1].data
    _mdp = grbdata['MDP 99%']
    t_obs = '100000'
    t_rep = '21600'
    plt.title('%i GRBs, $\Delta t_{obs}=%s s,$ $t_{repoint}=%s s$'\
              %(len(_mdp),t_obs,t_rep))
    plt.hist(_mdp*100, bins, alpha=0.5)
    plt.xlabel('2.-10. keV MDP (%)')
    plt.ylabel('Number of GRBs')
    overlay_tag()
    save_current_figure('all_grbs_MDP_histo', clear=False)

    # 2)----------------------------------------------------
    plt.figure(figsize=(10, 6), dpi=80)
    plt.title('%i GRBs, $\Delta t_{obs}=%s s,$ $t_{repoint}=%s s$'\
              %(len(_mdp),t_obs,t_rep))
    (n, bins, patches) = plt.hist(_mdp*100, bins, histtype='step', \
                                  cumulative=True)
    plt.xlabel('2.-10. keV MDP (%)')
    plt.ylabel('Cumulative number of GRBs')
    for i in range(0,30):
        print 'MDP %.2f%%: %i GRBs'%(i,n[i])
    overlay_tag()
    save_current_figure('all_grbs_MDP_cumulative', clear=False)

    # 3)------------------------------------------------------
    plt.figure(figsize=(10, 6), dpi=80)
    ax = plt.gca()
    _prompt_tstart = grbdata['PROMPT_START']
    _flux = grbdata['PROMPT_FLUX']
    _good_indexes = numpy.where(numpy.any(_prompt_tstart<350))
    print _good_indexes
    _flux = numpy.ma.masked_array(_flux, mask= _good_indexes)
    plt.scatter(_mdp*100, _flux, s=30, marker='.', color='gray')
    plt.xlabel('2.-10. keV MDP (%)')
    plt.ylabel('[erg $\cdot$ cm$^{-2}$]')
    plt.title('$\Delta t_{obs}=%s s,$ $t_{repoint}=%s s$'%(t_obs,t_rep))
    plt.xlim(1, 100)
    ax.set_yscale('log')
    ax.set_xscale('log')
    overlay_tag()
    save_current_figure('grb_MDP_prompt',clear=False)
    plt.show()
    """
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
    """

if __name__=='__main__':
    main()
