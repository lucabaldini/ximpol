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
from ximpol.config.crab_pulsar import pol_degree, pol_angle, pl_index,\
    pl_normalization


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'crab_pulsar.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'crab_pulsar')
SIM_DURATION = 100000.
NUM_PHASE_BINS = 20
E_BINNING = [1., 10.]

PIPELINE = xPipeline(clobber=False)


evt_file_path = '%s.fits' % OUT_FILE_PATH_BASE
analysis_file_path = '%s_analysis.txt' % OUT_FILE_PATH_BASE


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
    evt_file = xEventFile(evt_file_path)
    phase = evt_file.event_data['PHASE']
    bins = xEventBinningBase.equipopulated_binning(NUM_PHASE_BINS, phase, 0, 1)
    logger.info('Phase binning: %s' % bins)
    return bins


def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=evt_file_path)

# Since we need this in all the following steps, we create the binning at the
# top level of the module.
phasebins = _phase_binning()

def prepare():
    """Prepare the event data for the actual analysis.
    """
    for i, (_min, _max) in enumerate(zip(phasebins[:-1], phasebins[1:])):
        PIPELINE.xpselect(evt_file_path, phasemin=_min, phasemax=_max,
                          outfile=_sel_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='PHA1',
                       outfile=_pha1_file_path(i))

def analyze():
    """Analyze the data.
    """
    logger.info('Opening output file %s...' % analysis_file_path)
    analysis_file = open(analysis_file_path, 'w')
    for i, (_min, _max) in enumerate(zip(phasebins[:-1], phasebins[1:])):
        _mcube = xBinnedModulationCube(_mcube_file_path(i))
        _mcube.fit()
        _fit_results = _mcube.fit_results[0]
        _phase = 0.5*(_min + _max)
        _phase_err = 0.5*(_max - _min)
        _pol_deg = _fit_results.polarization_degree
        _pol_deg_err = _fit_results.polarization_degree_error
        _pol_angle = _fit_results.phase
        _pol_angle_err = _fit_results.phase_error
        _xspec_model = PIPELINE.xpxspec(_pha1_file_path(i))
        _index = _xspec_model(1).values[0]
        _index_err = _xspec_model(1).sigma
        # The division by the phase interval is a workaround and we should
        # keep track of that in xpselect.
        _norm = _xspec_model(2).values[0]/(_max - _min)
        _norm_err = _xspec_model(1).sigma/(_max - _min)
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
        _norm_err = numpy.loadtxt(analysis_file_path, unpack=True)
    plt.figure('Polarization degree')
    plt.errorbar(_phase, _pol_deg, xerr=_phase_err, yerr=_pol_deg_err, fmt='o',
                 label=sim_label)
    pol_degree.plot(show=False, label=mod_label)
    plt.axis([None, None, 0., 0.4])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('Polarization angle')
    plt.errorbar(_phase, _pol_angle, xerr=_phase_err, yerr=_pol_angle_err,
                 fmt='o', label=sim_label)
    pol_angle.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('PL normalization')
    plt.errorbar(_phase, _norm, xerr=_phase_err, yerr=_norm_err, fmt='o',
                 label=sim_label)
    pl_normalization.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.figure('PL index')
    plt.errorbar(_phase, _index, xerr=_phase_err, yerr=_index_err, fmt='o',
                 label=sim_label)
    pl_index.plot(show=False, label=mod_label)
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    plt.show()

def run():
    """Run all the tasks.
    """
    if os.path.exists(analysis_file_path):
        logger.info('%s exists, delete it if you want to recreate it.' %\
                    analysis_file_path)
    else:
        generate()
        prepare()
        analyze()
    plot()


if __name__ == '__main__':
    run()
