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

PHASE_BINS = [(0.,1.0), (0.03,0.33), (0.33, 0.43)]#, (0.95, 0.03)]

E_BINNING = [1.6, 6.]



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
    



def calcMDP():
    nebula_mcube_file = xBinnedModulationCube(NEBULA_MCUBE_FILE_PATH)
    nebula_counts = nebula_mcube_file.counts
    nebula_mdp =  nebula_mcube_file.mdp99
    txt = "Pulsar phase\t Emin - Emax\t Pulsar counts\t Nebula counts\t MDP\n"
    txt2 = "Pulsar phase\t Emin - Emax\t Pulsar MDP\t Nebula MDP\n"
    for i, (_min, _max) in enumerate(PHASE_BINS):
        pulse_diff = numpy.fabs(_max -_min)
        
        pulsar_phase_mcube_file = xBinnedModulationCube(_mcube_file_path(i))
        pulsar_phase_counts = pulsar_phase_mcube_file.counts
        pulsar_emean = pulsar_phase_mcube_file.emean
        pulsar_emin = pulsar_phase_mcube_file.emin
        pulsar_emax = pulsar_phase_mcube_file.emax
        pulsar_mdp = pulsar_phase_mcube_file.mdp99
        #scale the nebula counts for the time used for the pulsar phase
        scaled_nebula_counts = pulse_diff*nebula_counts
        
        count_sqrt = numpy.sqrt(pulsar_phase_counts + scaled_nebula_counts)
        #count_ratio = count_sqrt/pulsar_phase_counts
        eff_mu_pulsar =  pulsar_phase_mcube_file.effective_mu
        print "Phase:%s\tEff_mu_pulsar:%s\t pulsar counts:%s"%(PHASE_BINS[i],eff_mu_pulsar,pulsar_phase_counts)
   
        mdp = 4.292/eff_mu_pulsar*count_sqrt/pulsar_phase_counts
        
        
        txt += "%s\t %s - %s\t %s\t %s\t %.3f\n"%(PHASE_BINS[i], pulsar_emin[0], pulsar_emax[0], pulsar_phase_counts[0], scaled_nebula_counts[0], 100*mdp)

        txt2 +="%s\t %s - %s\t %.3f\t%.3f\n"%(PHASE_BINS[i], pulsar_emin[0], pulsar_emax[0],100*pulsar_mdp,100*nebula_mdp)
        
        
    print txt
    print
    print txt2
if __name__=='__main__':
    
    generate()
    prepare_pulsar()
    prepare_nebula()
    calcMDP()
