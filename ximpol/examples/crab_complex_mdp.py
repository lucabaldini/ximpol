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
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
NEBULA_FILE_PATH = '%s_nebula.fits' % OUT_FILE_PATH_BASE

SIM_DURATION = 10000

PHASE_BINS = [(0.03,0.33), (0.33, 0.43), (0.95, 0.03)]

E_BINNING = [1., 10.]
#OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def _sel_file_path(i):
    """Return the path to the i-th xpselct output file.
    """
    return '%s_phase%04d.fits' % (OUT_FILE_PATH_BASE, i)


def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)


def prepare_pulsar():
    """Prepare the event data for the actual analysis.
    """
    for i, (_min, _max) in enumerate(PHASE_BINS):
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,\
                    mcsrcid=1, rad=0.25, outfile=_sel_file_path(i))


def prepare_nebula():
    """Prepare the event data for the actual analysis.
    """
    PIPELINE.xpselect(EVT_FILE_PATH, mcsrcid=0, rad=0.25, \
                          outfile=NEBULA_FILE_PATH)


def calcMDP():
    nebula_file = xEventFile(NEBULA_FILE_PATH)
    nebula_counts = nebula_file.num_events()
    
    for i, (_min, _max) in enumerate(PHASE_BINS):
        pulse_diff = numpy.fabs(_max -_min)
        
        pulsar_phase_file = xEventFile(_sel_file_path(i))
        pulsar_phase_counts = pulsar_phase_file.num_events()
        scaled_nebula_counts = pulse_diff*nebula_counts
        
        count_sqrt = numpy.sqrt(pulsar_phase_counts + scaled_nebula_counts)
        count_ratio = count_sqrt/pulsar_phase_counts
        modf_ave = 0.3 #Just using a random value picked by hand
        mdp = (4.292/modf_ave)*count_ratio
        print "Pulsar phase: %s Pulsar counts:%s Nebula counts (in 15 arcsec) %s MDP:%s"%(PHASE_BINS[i], pulsar_phase_counts, scaled_nebula_counts, mdp)



if __name__=='__main__':
    
    generate()
    prepare_pulsar()
    prepare_nebula()
    calcMDP()
