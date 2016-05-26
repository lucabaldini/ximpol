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

CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_complex.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_complex')

EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE

REG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'crab_nebula_9regions.reg')
MAP_FILE_PATH = os.path.join(XIMPOL_DATA, 'crab_complex_cmap.fits')


DURATION = 10000
E_BINNING = [1., 4., 6.]

pipeline = xPipeline(clobber=False)

def get_sel_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'crab_complex_reg%04d.fits' % i)

def get_mcube_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'crab_complex_reg%04d_mcube.fits' % i)

def generate():
    """
    """
    pipeline.xpobssim(configfile=CFG_FILE_PATH, duration=DURATION,
                      outfile=EVT_FILE_PATH)

def select_and_bin():
    """
    """
    logger.info('Creating the mapcube for the entire source...')
    pipeline.xpbin(EVT_FILE_PATH, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING)
    logger.info('Opening region file %s...' % REG_FILE_PATH)
    regions = pyregion.open(REG_FILE_PATH)
    logger.info('Found %d regions...' % len(regions))
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        rad *= 60.
        logger.info('Analyzing region at ra = %s, dec = %s' % (ra, dec))
        sel_file_path = get_sel_file_path(i)
        mcube_file_path = get_mcube_file_path(i)
        pipeline.xpselect(EVT_FILE_PATH, ra=ra, dec=dec, rad=rad,
                          outfile=sel_file_path)
        pipeline.xpbin(sel_file_path, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile = mcube_file_path)

def plot(save=False):
    logger.info('Plotting stuff...')
    pipeline.xpbin(EVT_FILE_PATH, algorithm='CMAP', outfile=MAP_FILE_PATH)
    regions = pyregion.open(REG_FILE_PATH)
    full_map = xBinnedMap(MAP_FILE_PATH)
    fig_all = full_map.plot(show=False)

    #First make the global plot
    fig_global = full_map.plot(show=False, subplot=(1, 2, 1))
    plt.subplots_adjust(hspace=0.001)
    global_mcube_file_path = os.path.join(XIMPOL_DATA, 'crab_complex_mcube.fits' )
    mcube = xBinnedModulationCube(global_mcube_file_path)
    mcube.plot(show=False, analyze=False, xsubplot=1)
    fig_global.save(os.path.join(XIMPOL_DATA, 'crab_complex_global.png'))
    plt.clf()
    
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        #fig_all.show_circles(ra, dec, rad, lw=1)

        fig = full_map.plot(show=False, subplot=(1, 2, 1))
        plt.subplots_adjust(hspace=0.001)
        fig.show_circles(ra, dec, rad, lw=1)
        mcube_file_path = get_mcube_file_path(i)
        mcube = xBinnedModulationCube(mcube_file_path)
        mcube.plot(show=False, analyze=False, xsubplot=1)

        scale_x  = rad/numpy.cos(numpy.deg2rad(dec)) # This is to take into account the effect of the projection.
        scale_y  = rad

        for j,fit in enumerate(mcube.fit_results):
            angle = fit.phase
            angle_error = fit.phase_error
            degree = fit.polarization_degree

            #nangles=20
            #for t in range(nangles):
            #    dx = scale_x*numpy.cos(numpy.pi*2.0*t/float(nangles))#angle)
            #    dy = scale_y*numpy.sin(numpy.pi*2.0*t/float(nangles))#angle)
            #    fig.show_arrows(ra, dec, dx, dy, color='w', alpha=1, width=1,head_width=0, head_length=0)
            #    pass

            dx = scale_x*numpy.cos(angle)
            dy = scale_y*numpy.sin(angle)

            dx1 = scale_x*degree*numpy.cos(angle+angle_error)
            dy1 = scale_y*degree*numpy.sin(angle+angle_error)
            dx2 = scale_x*degree*numpy.cos(angle-angle_error)
            dy2 = scale_y*degree*numpy.sin(angle-angle_error)

            fig.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx, -dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, dx1, dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, dx2, dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)

            fig_all.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx, -dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, dx1, dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, dx2, dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)


        if save:
            fig.save(mcube_file_path.replace('.fits', '.png'))
            plt.clf()
        else:
            plt.show()
        pass
    fig_all.save(os.path.join(XIMPOL_DATA, 'crab_complex_reg_all.png'))
    
    
def plot_doc():
    """
    """
    img = xFITSImage(os.path.join(XIMPOL_DATA, 'crab_complex_cmap.fits'))
    fig = img.plot(show=False)
    xFITSImage.add_label(fig, 'XIPE %d ks' %DURATION/1000.)
    plt.show()


if __name__ == '__main__':
    generate()
    select_and_bin()
    plot(True)
    #plot_doc()




import os
import numpy

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedModulationCube, xEventBinningBase
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.config.crab_pulsar import pol_degree_spline, pol_angle_spline,\
    pl_index_spline, pl_normalization_spline


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_complex.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_complex')

#PULSAR_CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_pulsar.py')
#PULSAR_OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_pulsar')

#NEBULA_CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_nebula.py')
#NEBULA_OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_nebula')

EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE

#NEBULA_EVT_FILE_PATH = '%s.fits' % NEBULA_OUT_FILE_PATH_BASE

NEBULA_SELECTED_FILE_PATH = '%s_nebula_selected.fits' % OUT_FILE_PATH_BASE
NEBULA_MCUBE_FILE_PATH = '%s_nebula_mcube.fits' % OUT_FILE_PATH_BASE

SIM_DURATION = 10000

PHASE_BINS = [(0.,0.05), (0.05,0.25), (0.25, 0.45), (0.45, 0.9), (0.9,1.0)]
#NUM_PHASE_BINS = 25
#PHASE_BINS = numpy.linspace(0., 1., NUM_PHASE_BINS)

E_BINNING = [1.0, 3., 5., 8., 10.]
#E_BINNING = numpy.linspace(1., 10., 5)

OUTPUT_FOLDER = '/data/work/ximpol/ximpol/examples/crab/'

MDP_OUTPUT_FILE = os.path.join(OUTPUT_FOLDER,'MDP_CrabPulsar_imaging_fE.txt')

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


def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)

   

def prepare_pulsar():
    """Prepare the event data for the actual analysis.
    """
    for i, (_min, _max) in enumerate(PHASE_BINS):
            #zip(PHASE_BINS[:-1],PHASE_BINS[1:])):
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min,
                          phasemax=_max, mcsrcid=1,rad=0.25,
                          outfile=_sel_file_path(i))
        
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))

def prepare_nebula():
    """Prepare the event data for the actual analysis.
    """
    PIPELINE.xpselect(EVT_FILE_PATH, mcsrcid=0, rad=0.25, \
                          outfile=NEBULA_SELECTED_FILE_PATH)
    
    PIPELINE.xpbin(NEBULA_SELECTED_FILE_PATH, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=NEBULA_MCUBE_FILE_PATH)
    



def calcMDP(plot=False):
    mdp_file = open(MDP_OUTPUT_FILE,'w')
    mdp_file.write('#Simulation time %s sec \n'%SIM_DURATION)
    mdp_file.write('#Phase Ave delta phase Mean energy  delta energy mdp99\n')
    nebula_mcube_file = xBinnedModulationCube(NEBULA_MCUBE_FILE_PATH)
    nebula_counts = nebula_mcube_file.counts
    nebula_mdp =  nebula_mcube_file.mdp99
    txt = "Pulsar phase\t Emin - Emax\t Pulsar counts\t Nebula counts\t MDP\n"
    MDP99_PULSAR = []
    phase = []
    phase_err = [] 
    for i, (_min, _max) in enumerate(PHASE_BINS):
            #zip(PHASE_BINS[:-1],PHASE_BINS[1:])):

        pulse_diff = numpy.fabs(_max -_min)
        _phase_ave = 0.5*(_min + _max)
        phase.append(_phase_ave)
        _phase_err = 0.5*(_max - _min)
        phase_err.append(_phase_err)
        
        pulsar_phase_mcube_file = xBinnedModulationCube(_mcube_file_path(i))
        
        
        pulsar_emean = pulsar_phase_mcube_file.emean
        for j, _energy_mean in enumerate(pulsar_emean):
            pulsar_emin = pulsar_phase_mcube_file.emin[j]
            pulsar_emax = pulsar_phase_mcube_file.emax[j]
            pulsar_e_err = 0.5*(pulsar_emax-pulsar_emin)
            
            pulsar_phase_counts = pulsar_phase_mcube_file.counts[j]
            pulsar_mdp = pulsar_phase_mcube_file.mdp99[j]
            #scale the nebula counts for the time used for the pulsar phase
            scaled_nebula_counts = pulse_diff*nebula_counts[j]
        
            count_sqrt = numpy.sqrt(pulsar_phase_counts + scaled_nebula_counts)
       
            eff_mu_pulsar =  pulsar_phase_mcube_file.effective_mu[j]
            
            mdp = 4.292/eff_mu_pulsar*count_sqrt/pulsar_phase_counts
            MDP99_PULSAR.append(100*mdp)
            _data = (_phase_ave, _phase_err, _energy_mean, pulsar_e_err, pulsar_e_err, mdp)
            
            _fmt = ('%.4e   ' * len(_data)).strip()
            _fmt = '%s\n' % _fmt
            _line = _fmt % _data
            mdp_file.write(_line)
           
            #txt += "%s\t %s - %s\t %s\t %s\t %.3f\n"%(PHASE_BINS[i], pulsar_emin[0], pulsar_emax[0], pulsar_phase_counts[0], scaled_nebula_counts[0], 100*mdp)
        
    MDP99_PULSAR = numpy.array(MDP99_PULSAR)
    PHASE = numpy.array(phase)
    mdp_file.close()
    if plot:
        scale_factor = 10
        sim_label = 'XIPE %s ks' % (SIM_DURATION*scale_factor/1000.)
        lc_label = 'Light curve'
        plt.errorbar(PHASE, MDP99_PULSAR*(1/numpy.sqrt(10)),xerr=phase_err, label=sim_label,fmt='o')
        pl_normalization_spline.plot(scale=10., show=False,
                                color='lightgray',label=lc_label)
        plt.ylabel('MDP 99\%')
        plt.legend()
        plt.savefig('crab_complex_mdp_nonimaging_%i.png'%(SIM_DURATION*scale_factor/1000.))
      
        #plt.show()
    print txt
   


def makeMDPComparisonPlot(imaging_file_path):
    
    non_imaging_file_path = imaging_file_path.replace('_imaging.txt','_non_imaging.txt')
    print "Using %s for non imaging file path"%non_imaging_file_path
    
    _phase_ave, _phase_err, mdp_imaging = numpy.loadtxt(imaging_file_path, unpack=True)
    _phase_ave, _phase_err, mdp_nonimaging = numpy.loadtxt(non_imaging_file_path, unpack=True)
    scale_factor = 10.
    print "Improvement with imaging"
    for i, phase in enumerate(_phase_ave):
        imaging = 100*mdp_imaging[i]*(1/numpy.sqrt(scale_factor))
        non_imaging = 100*mdp_nonimaging[i]*(1/numpy.sqrt(scale_factor))
        print "%s\t Imaging (non):%s (%s) %s"%(phase,imaging,non_imaging,non_imaging/imaging)
        
   
    sim_label_imaging = 'XIPE %s ks\n Imaging 15"' % (SIM_DURATION*scale_factor/1000.)
    sim_label_nonimaging = 'Non imaging'
    lc_label = 'Light curve'
    plt.errorbar(_phase_ave, 100*mdp_imaging*(1/numpy.sqrt(scale_factor)),xerr=_phase_err, label=sim_label_imaging,fmt='o',markersize=6)
    
    plt.errorbar(_phase_ave, 100*mdp_nonimaging*(1/numpy.sqrt(scale_factor)),xerr=_phase_err, label=sim_label_nonimaging,fmt='v',markersize=6)
    
    pl_normalization_spline.plot(scale=10., show=False,
                                color='darkgray',label=lc_label)
    #on_phase = 0.25, 0.45
    #off_phase = 0.45,0.9
    plt.axvspan(0.25, 0.45, color='r', alpha=0.45, lw=0)
    plt.axvspan(0.45, 0.9, color='gray', alpha=0.25, lw=0)
    
    plt.ylabel('MDP 99\%')
    plt.legend()
    plt.savefig('crab_complex_mdp_imaging_vs_nonimaging_%i_shaded.png'%(SIM_DURATION*scale_factor/1000.))
    plt.show()


def makeMDP_fE_ComparisonPlot(file_path):
    scale_factor = 10.
    
    (_phase_ave, _phase_err, _energy_mean, pulsar_e_err, pulsar_e_err, mdp) =\
                                                numpy.loadtxt(file_path, unpack=True)
    print "Phase ave:",_phase_ave
    print
    print "Energy mean", _energy_mean
    #phase_values = [0.025, 0.15, 0.35, 0.675, 0.95]
    phase_values = [0.35,0.675]
    on, on_phase_color = (0.35,'r')
    off, off_phase_color = (0.675,'gray') 
    #for phase in phase_values:
        
    plt.errorbar(_energy_mean[_phase_ave==on], 100*mdp[_phase_ave==on]*(1/numpy.sqrt(scale_factor)),xerr=pulsar_e_err[_phase_ave==on], label='On Phase',fmt='o',markersize=6,ls='--',color=on_phase_color)

    plt.errorbar(_energy_mean[_phase_ave==off], 100*mdp[_phase_ave==off]*(1/numpy.sqrt(scale_factor)),xerr=pulsar_e_err[_phase_ave==off], label='Off Phase',fmt='o',markersize=6,ls='--',color=off_phase_color)
    
    plt.legend()
    plt.ylabel('MPD 99\%')
    plt.xlabel('Energy (keV)')
    plt.savefig('crab_complex_mdp_imaging_fE_%i.png'%(SIM_DURATION*scale_factor/1000.))
    plt.show()
    
    
if __name__=='__main__':
    
    generate()
    #prepare_pulsar()
    #prepare_nebula()
    #calcMDP(plot=False)
    #makeMDPComparisonPlot('MDP_CrabPulsar_imaging.txt')
    #makeMDP_fE_ComparisonPlot('MDP_CrabPulsar_imaging_fE.txt')
