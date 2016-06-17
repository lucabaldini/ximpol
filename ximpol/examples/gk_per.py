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
from ximpol.config.gk_per import pol_degree_spline, phasogram_spline,\
    GK_PER_POLARIZATION_ANGLE


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'gk_per.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'gk_per')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
SIM_DURATION = 800000.
PHASE_BINNING = numpy.linspace(0., 1., 9)
PHASE_BINNING = zip(PHASE_BINNING[:-1], PHASE_BINNING[1:])
E_BINNING = [2., 8.]


def _sel_file_path(i):
    """Return the path to the i-th xpselect output file.
    """
    return '%s_phase%04d.fits' % (OUT_FILE_PATH_BASE, i)

def _mcube_file_path(i):
    """Return the path to the i-th xpbin MCUBE output file.
    """
    return '%s_phase%04d_mcube.fits' % (OUT_FILE_PATH_BASE, i)


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)

def bin():
    """Bin the event file in different flavors for the actual analysis.
    """
    PIPELINE.xpbin(EVT_FILE_PATH, algorithm='CMAP')
    for i, (_min, _max) in enumerate(PHASE_BINNING):
        PIPELINE.xpselect(EVT_FILE_PATH, phasemin=_min, phasemax=_max,
                          outfile=_sel_file_path(i))
        PIPELINE.xpbin(_sel_file_path(i), algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile=_mcube_file_path(i))

def analyze():
    """
    """
    if os.path.exists(ANALYSIS_FILE_PATH):
        logger.info('%s exists, delete it if you want to recreate it.' %\
                    ANALYSIS_FILE_PATH)
        return
    logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
    for i, (_min, _max) in enumerate(PHASE_BINNING):
        _mcube = xBinnedModulationCube(_mcube_file_path(i))
        _mcube.fit()
        _fit_results = _mcube.fit_results[0]
        print _fit_results
        _phase = 0.5*(_min + _max)
        _phase_err = 0.5*(_max - _min)
        _pol_deg = _fit_results.polarization_degree
        _pol_deg_err = _fit_results.polarization_degree_error
        _pol_angle = _fit_results.phase
        _pol_angle_err = _fit_results.phase_error
        _data = (_phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,
                 _pol_angle_err)
        _fmt = ('%.4e   ' * len(_data)).strip()
        _fmt = '%s\n' % _fmt
        _line = _fmt % _data
        analysis_file.write(_line)
    analysis_file.close()

def plot(save_plots=False):
    """
    """
    sim_label = 'XIPE %s ks' % (SIM_DURATION/1000.)
    mod_label = 'Input model'
    _phase, _phase_err, _pol_deg, _pol_deg_err, _pol_angle,\
        _pol_angle_err = numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    _colors = ['blue']*len(_pol_deg)
    plt.figure('Polarization degree')
    _good_fit = _pol_deg > 2*_pol_deg_err
    _bad_fit = numpy.logical_not(_good_fit)
    plt.errorbar(_phase[_good_fit], _pol_deg[_good_fit],
                 xerr=_phase_err[_good_fit], yerr=_pol_deg_err[_good_fit],
                 fmt='o', label=sim_label, color='blue')
    plt.errorbar(_phase[_bad_fit], _pol_deg[_bad_fit],
                 xerr=_phase_err[_bad_fit], yerr=_pol_deg_err[_bad_fit],
                 fmt='o', label=sim_label, color='gray')
    pol_degree_spline.plot(show=False, label=mod_label, color='green')
    plt.axis([0., 1., 0., 0.1])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save_plots:
        plt.savefig('gk_per_polarization_degree.png')
    plt.figure('Polarization angle')
    plt.errorbar(_phase[_good_fit], _pol_angle[_good_fit],
                 xerr=_phase_err[_good_fit], yerr=_pol_angle_err[_good_fit],
                 fmt='o', label=sim_label, color='blue')
    plt.errorbar(_phase[_bad_fit], _pol_angle[_bad_fit],
                 xerr=_phase_err[_bad_fit], yerr=_pol_angle_err[_bad_fit],
                 fmt='o', label=sim_label, color='gray')
    _x = numpy.linspace(0., 1., 100)
    _y = numpy.full(100, GK_PER_POLARIZATION_ANGLE)
    plt.plot(_x, _y, label=mod_label, color='green')
    plt.axis([0., 1., 0., 1.5])
    plt.xlabel('Rotational phase')
    plt.ylabel('Polarization angle [rad]')
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save_plots:
        plt.savefig('gk_per_polarization_angle.png')
    plt.show()

def run(save_plots):
    """
    """
    generate()
    bin()
    analyze()
    plot(save_plots)


if __name__ == '__main__':
    run(save_plots=True)
