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

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedModulationCube, xEventBinningBase
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.config.crab_pulsar import pol_degree_spline, pol_angle_spline,\
    pl_index_spline, pl_normalization_spline


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_pulsar.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_pulsar')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
SIM_DURATION = 100000.
NUM_PHASE_BINS = 20
PHASE_BINNING = None
E_BINNING = [1., 10.]


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

def _phase_binning():
    """Read the input event file and create an equipopulated binning in the
    pulsar phase.
    """
    evt_file = xEventFile(EVT_FILE_PATH)
    phase = evt_file.event_data['PHASE']
    bins = xEventBinningBase.equipopulated_binning(NUM_PHASE_BINS, phase, 0, 1)
    logger.info('Phase binning: %s' % bins)
    return bins


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
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,
                          outfile=_sel_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='PHA1',
                       outfile=_pha1_file_path(i))

def analyze():
    """Analyze the data.
    """
    logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
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
        _spec_fitter = PIPELINE.xpxspec(_pha1_file_path(i), plot=False)
        (_index, _index_err), (_norm, _norm_err) =_spec_fitter.fit_parameters()
        # The division by the phase interval is a workaround and we should
        # keep track of that in xpselect.
        _norm /= (_max - _min)
        _norm_err /= (_max - _min)
        _data = (_phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,
                 _pol_angle_err, _index, _index_err, _norm, _norm_err)
        _fmt = ('%.4e   ' * len(_data)).strip()
        _fmt = '%s\n' % _fmt
        _line = _fmt % _data
        analysis_file.write(_line)
    analysis_file.close()

def plot():
    """Plot the stuff in the analysis file.
    """
    sim_label = 'XIPE %s ks' % (SIM_DURATION/1000.)
    mod_label = 'Input model'
    _phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,\
        _pol_angle_err, _index, _index_err, _norm,\
        _norm_err = numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    plt.figure('Polarization degree')
    plt.errorbar(_phase, _pol_deg, xerr=_phase_err, yerr=_pol_deg_err, fmt='o',
                 label=sim_label)
    pol_degree_spline.plot(show=False, label=mod_label)
    plt.axis([None, None, 0., 0.4])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('Polarization angle')
    plt.errorbar(_phase, _pol_angle, xerr=_phase_err, yerr=_pol_angle_err,
                 fmt='o', label=sim_label)
    pol_angle_spline.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('PL normalization')
    plt.errorbar(_phase, _norm, xerr=_phase_err, yerr=_norm_err, fmt='o',
                 label=sim_label)
    pl_normalization_spline.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('PL index')
    plt.errorbar(_phase, _index, xerr=_phase_err, yerr=_index_err, fmt='o',
                 label=sim_label)
    pl_index_spline.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.show()

def run():
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
    plot()


if __name__ == '__main__':
    run()
