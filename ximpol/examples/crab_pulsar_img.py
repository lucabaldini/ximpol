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
import numpy
from astropy.io import fits
from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol import xpColor

from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedModulationCube, xEventBinningBase, xBinnedMap
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.config.crab_pulsar import pol_degree_spline, pol_angle_spline,\
    pl_index_spline, pl_normalization_spline


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_pulsar.py')
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_complex.py')
#CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_nebula.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_complex_zero')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
SIM_DURATION = 100000.
NUM_PHASE_BINS = 25
EQP_BINNING = False
PHASE_BINNING = None
E_BINNING = [1., 10.]
OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')
RA=83.633083
DEC=22.0145
RAD_MAP=2.0 #arcmin
RAD_ANA=15.0/60. #arcmin

"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def _sel_file_path(i):
    """Return the path to the i-th xpselct output file.
    """
    return '%s_phase%04d.fits' % (OUT_FILE_PATH_BASE, i)

def _mcube_file_path(i):
    """Return the path to the i-th xpbin MCUBE output file.
    """
    return '%s_phase%04d_mcube.fits' % (OUT_FILE_PATH_BASE, i)

def _pha1_file_path(i):
    """Return the path to the i-th xpbin PHA1 output file.
    """
    return '%s_phase%04d_pha1.fits' % (OUT_FILE_PATH_BASE, i)

def _cmap_file_path(i):
    """Return the path to the i-th xpbin PHA1 output file.
    """
    return '%s_phase%04d_cmap.fits' % (OUT_FILE_PATH_BASE, i)

def _phase_binning():
    """Read the input event file and create an equipopulated binning in the
    pulsar phase.
    """
    if EQP_BINNING:
        evt_file = xEventFile(EVT_FILE_PATH)
        phase = evt_file.event_data['PHASE']
        return xEventBinningBase.equipopulated_binning(NUM_PHASE_BINS, phase,
                                                       0., 1.)
    else:
        return numpy.linspace(0., 1., NUM_PHASE_BINS)
    
def _plot_arrows(ra,dec,rad,angle,angle_error,degree,fig,color='r'):
    scale_x  = rad/numpy.cos(numpy.deg2rad(dec)) # This is to take into account the effect of the projection.
    scale_y  = rad
    
    dx = scale_x*numpy.cos(angle)
    dy = scale_y*numpy.sin(angle)
    
    dx1 = scale_x*degree*numpy.cos(angle+angle_error)
    dy1 = scale_y*degree*numpy.sin(angle+angle_error)
    dx2 = scale_x*degree*numpy.cos(angle-angle_error)
    dy2 = scale_y*degree*numpy.sin(angle-angle_error)
    
    fig.show_arrows(ra, dec, dx, dy, color=color, alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
    fig.show_arrows(ra, dec, -dx, -dy, color=color, alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
    fig.show_arrows(ra, dec, dx1, dy1, color=color, alpha=1, width=1,head_width=0, head_length=0)
    fig.show_arrows(ra, dec, -dx1, -dy1, color=color, alpha=1, width=1,head_width=0, head_length=0)
    fig.show_arrows(ra, dec, dx2, dy2, color=color, alpha=1, width=1,head_width=0, head_length=0)
    fig.show_arrows(ra, dec, -dx2, -dy2, color=color, alpha=1, width=1,head_width=0, head_length=0)
    pass
    
def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)

def prepare():
    """Prepare the event data for the actual analysis.
    """
    for i, (_min, _max) in enumerate(zip(PHASE_BINNING[:-1],
                                         PHASE_BINNING[1:])):

        #1) select on the phse & large field of view:
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,
                          ra=RA, dec=DEC, rad=RAD_MAP,
                          outfile=_sel_file_path(i))
        #2) make a CMAP
        PIPELINE.xpbin(_sel_file_path(i), algorithm='CMAP',
                       outfile=_cmap_file_path(i))
        #3) select on the ROI:
        os.system('rm %s' % _sel_file_path(i)) 
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,
                          ra=RA, dec=DEC, rad=RAD_ANA,
                          outfile=_sel_file_path(i))
        #4) make another small cmap
        PIPELINE.xpbin(_sel_file_path(i), algorithm='CMAP',
                       outfile=_cmap_file_path(i).replace('_cmap.fits','_cmap_sel.fits'))

        
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='PHA1',
                       outfile=_pha1_file_path(i))

def analyze():
    """Analyze the data.
    """
    logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
    _nevents=[]
    for i, (_min, _max) in enumerate(zip(PHASE_BINNING[:-1],
                                         PHASE_BINNING[1:])):
        _nevents.append(fits.open(_sel_file_path(i))['EVENTS'].header['NAXIS2'])
        logger.info("PHASE BIN %d [%.1f-%.1f] - EVENTS = %d" %(i,_min,_max,_nevents[i]))        
        pass
    out_phase=numpy.min(_nevents)
    for i, (_min, _max) in enumerate(zip(PHASE_BINNING[:-1],
                                         PHASE_BINNING[1:])):
        _mcube = xBinnedModulationCube(_mcube_file_path(i))
        _mcube.fit()
        _fit_results = _mcube.fit_results[0]
        _phase = 0.5*(_min + _max)
        _phase_err = 0.5*(_max - _min)
        _pol_deg = _fit_results.polarization_degree
        _pol_deg_err = _fit_results.polarization_degree_error
        _pol_angle = _fit_results.phase
        _pol_angle_err = _fit_results.phase_error
        if (_nevents[i]-out_phase)>0:
            _factor=(_nevents[i])/(_nevents[i]-out_phase)
            print '------->',_nevents[i],out_phase,_nevents[i]-out_phase,_factor
            _pol_deg*=_factor
            _pol_deg_err*=_factor
        else:
            _pol_deg=1e-6
            _pol_deg_err=1e-6
            pass
        if _pol_deg>1.0: _pol_deg=1.0
        #_total_px=_pol_deg*numpy.cos(_pol_angle)
        #_total_py=_pol_deg*numpy.sin(_pol_angle)
        #_nebula_px=0.157*numpy.cos(numpy.radians(161.1))
        #_nebula_py=0.157*numpy.sin(numpy.radians(161.1))
        #_pulsar_px=_total_px-_nebula_px
        #_pulsar_py=_total_py-_nebula_py
        #_pol_angle=numpy.arctan2(_pulsar_py,_pulsar_px)
        #_pol_deg=numpy.sqrt(_pulsar_px*_pulsar_px+_pulsar_py*_pulsar_py)
        
        if False:
            _spec_fitter = PIPELINE.xpxspec(_pha1_file_path(i), plot=False)
            (_index, _index_err), (_norm, _norm_err) = _spec_fitter.fit_parameters()
            # The division by the phase interval is a workaround and we should
            # keep track of that in xpselect.
            _norm /= (_max - _min)
            _norm_err /= (_max - _min)
        else:
            _norm=0.0
            _norm_err=0.0
            _index=0.0
            _index_err=0.0
            pass
        
        _data = (_phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,
                 _pol_angle_err, _index, _index_err, _norm, _norm_err)
        _fmt = ('%.4e   ' * len(_data)).strip()
        _fmt = '%s\n' % _fmt
        _line = _fmt % _data
        analysis_file.write(_line)
        ### Plot the cmap with the arrow...
        full_map = xBinnedMap(_cmap_file_path(i))
        full_map.image.vmin=0
        full_map.image.vmax=10000        
        fig_map = full_map.plot(show=False)
        fig_map.show_circles(RA, DEC, RAD_ANA/60.0, lw=1)
                
        _plot_arrows(ra=RA,dec=DEC,rad=RAD_ANA/60.,
                     angle=_phase,angle_error=_phase_err,
                     degree=_pol_deg,fig=fig_map,color='b')
        
        fig_map_file=_cmap_file_path(i).replace('.fits','.png')        
        print 'saving map in...',fig_map_file
        fig_map.save(fig_map_file)
        ###
    analysis_file.close()

def plot(save=False):
    """Plot the stuff in the analysis file.
    """
    sim_label = 'XIPE %s ks' % (SIM_DURATION/1000.)
    sim_label = 'OBS: %s ks' % (SIM_DURATION/1000.)
    mod_label = 'Input model'
    lc_label = 'Light curve'
    _phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,\
        _pol_angle_err, _index, _index_err, _norm,\
        _norm_err = numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    plt.figure('Polarization degree')
    pl_normalization_spline.plot(scale=0.12, show=False, color='lightgray',
                                 label=lc_label)

    _pol_deg_err_p=numpy.where((_pol_deg+_pol_deg_err)<1.0,_pol_deg_err,1.0-_pol_deg)
    _pol_deg_err_m=numpy.where((_pol_deg-_pol_deg_err)>0.0,_pol_deg_err,_pol_deg)
    #if (_pol_deg+_pol_deg_err)>1.0: _pol_deg_err=1.0-yerr_p
    
    plt.errorbar(_phase, _pol_deg, xerr=_phase_err, yerr=[_pol_deg_err_m,_pol_deg_err_p], fmt='o',
                 label=sim_label)
    pol_degree_spline.plot(show=False, label=mod_label)
    plt.axis([0., 1., 0., 0.5])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('crab_polarization_degree', OUTPUT_FOLDER, False)
    plt.figure('Polarization angle')
    pl_normalization_spline.plot(scale=0.4, offset=1.25, show=False,
                                 color='lightgray', label=lc_label)
    plt.errorbar(_phase, _pol_angle, xerr=_phase_err, yerr=_pol_angle_err,
                 fmt='o', label=sim_label)
    pol_angle_spline.plot(show=False, label=mod_label)
    plt.axis([0., 1., 1.25, 3.])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('crab_polarization_angle', OUTPUT_FOLDER, False)
    plt.figure('PL normalization')
    plt.errorbar(_phase, _norm, xerr=_phase_err, yerr=_norm_err, fmt='o',
                 label=sim_label)
    pl_normalization_spline.plot(show=False, label=mod_label)
    plt.axis([0., 1., None, None])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('crab_pl_norm', OUTPUT_FOLDER, False)
    plt.figure('PL index')
    pl_normalization_spline.plot(scale=0.18, offset=1.3, show=False,
                                 color='lightgray', label=lc_label)
    plt.errorbar(_phase, _index, xerr=_phase_err, yerr=_index_err, fmt='o',
                 label=sim_label)
    pl_index_spline.plot(show=False, label=mod_label)
    plt.axis([0., 1., 1.3, 2.1])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('crab_pl_index', OUTPUT_FOLDER, False)
    plt.show()

def run(save_plots=False):
    """Run all the tasks.
    """
    if os.path.exists(ANALYSIS_FILE_PATH):
        logger.info('%s exists, delete it if you want to recreate it.' %\
                    ANALYSIS_FILE_PATH)
    else:
        generate()
        global PHASE_BINNING
        PHASE_BINNING = _phase_binning()
        prepare()
        analyze()
    plot(save_plots)


if __name__ == '__main__':
    run(save_plots=True)
