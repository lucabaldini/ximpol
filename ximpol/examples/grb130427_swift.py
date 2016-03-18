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

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_EXAMPLES, XIMPOL_DOC
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedMap, xBinnedModulationCube
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.config.grb130427_swift import pol_degree_spline,\
    polarization_angle, PL_INDEX
from ximpol.utils.matplotlib_ import save_current_figure
from ximpol.evt.event import xEventFile
from ximpol.core.spline import xInterpolatedUnivariateSpline


"""Script-wide simulation and analysis settings.
"""
CFG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'grb130427_swift.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'grb130427_swift')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
SIM_DURATION = 1000000.
NUM_TIME_BINS = 25
TIME_BINNING = numpy.logspace(numpy.log10(131.2489717),
                              numpy.log10(1000131.2489717), NUM_TIME_BINS)
E_BINNING = [1., 10.]
OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def _sel_file_path(i):
    """Return the path to the i-th xpselct output file.
    """
    return '%s_interval%04d.fits' % (OUT_FILE_PATH_BASE, i)

def _mcube_file_path(i):
    """Return the path to the i-th xpbin MCUBE output file.
    """
    return '%s_interval%04d_mcube.fits' % (OUT_FILE_PATH_BASE, i)

def _pha1_file_path(i):
    """Return the path to the i-th xpbin PHA1 output file.
    """
    return '%s_interval%04d_pha1.fits' % (OUT_FILE_PATH_BASE, i)


def generate():
    """Generate the events.
    """
    PIPELINE.xpobssim(configfile=CFG_FILE_PATH, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)

def prepare():
    """Prepare the event data for the actual analysis.
    """
    PIPELINE.xpbin(EVT_FILE_PATH, algorithm='LC', tbinalg='LOG')
    for i, (_min, _max) in enumerate(zip(TIME_BINNING[:-1],
                                         TIME_BINNING[1:])):
        PIPELINE.xpselect(EVT_FILE_PATH, tmin=_min, tmax=_max,
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
    for i, (_min, _max) in enumerate(zip(TIME_BINNING[:-1],
                                         TIME_BINNING[1:])):
        _mcube = xBinnedModulationCube(_mcube_file_path(i))
        _mcube.fit()
        _fit_results = _mcube.fit_results[0]
        _sel_file = xEventFile(_sel_file_path(i))
        _time = numpy.average(_sel_file.event_data['TIME'])
        _sel_file.close()
        _time_errp = _max - _time
        _time_errm = _time - _min
        _pol_deg = _fit_results.polarization_degree
        _pol_deg_err = _fit_results.polarization_degree_error
        _pol_angle = _fit_results.phase
        _pol_angle_err = _fit_results.phase_error
        _spec_fitter = PIPELINE.xpxspec(_pha1_file_path(i), plot=False)
        (_index, _index_err), (_norm, _norm_err) = _spec_fitter.fit_parameters()
        # The division by the phase interval is a workaround and we should
        # keep track of that in xpselect.
        _norm /= (_max - _min)
        _norm_err /= (_max - _min)
        _data = (_time, _time_errp, _time_errm, _pol_deg, _pol_deg_err,
                 _pol_angle, _pol_angle_err, _index, _index_err, _norm,
                 _norm_err)
        _fmt = ('%.4e   ' * len(_data)).strip()
        _fmt = '%s\n' % _fmt
        _line = _fmt % _data
        analysis_file.write(_line)
    analysis_file.close()

def plot(save=False):
    """Plot the stuff in the analysis file.
    """
    sim_label = 'XIPE'
    mod_label = 'Input model'
    lc_label = 'Light curve'
    _time, _time_errp, _time_errm, _pol_deg, _pol_deg_err, _pol_angle,\
        _pol_angle_err, _index, _index_err, _norm,\
        _norm_err = numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    _pol_angle = numpy.degrees(_pol_angle)
    _pol_angle_err = numpy.degrees(_pol_angle_err)
    plt.figure('Polarization degree')
    pol_degree_spline.plot(show=False, label=mod_label, logx=True)
    plt.errorbar(_time, _pol_deg, xerr=[_time_errm, _time_errp],
                 yerr=_pol_deg_err, fmt='o', label=sim_label)
    plt.axis([100., 1e6, 0., 0.6])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('grb130427_swift_polarization_degree',
                            OUTPUT_FOLDER, False)
    plt.figure('Polarization angle')
    _x = pol_degree_spline.x
    _y = numpy.full(len(_x), polarization_angle(_x, None, None, None))
    _y = numpy.degrees(_y)
    fmt = dict(xname='Time', xunits='s', yname='Polarization angle',
               yunits=r'$^\circ$')
    _s = xInterpolatedUnivariateSpline(_x, _y, **fmt)
    _s.plot(logx=True, show=False, label=mod_label)
    plt.errorbar(_time, _pol_angle, xerr=[_time_errm, _time_errp],
                 yerr=_pol_angle_err, fmt='o', label=sim_label)
    plt.axis([100., 1e6, None, None])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('grb130427_swift_polarization_angle',
                            OUTPUT_FOLDER, False)
    plt.figure('PL index')
    _y = numpy.full(len(_x), PL_INDEX)
    fmt = dict(xname='Time', xunits='s', yname='PL index')
    _s = xInterpolatedUnivariateSpline(_x, _y, **fmt)
    _s.plot(logx=True, show=False, label=mod_label)
    plt.errorbar(_time, _index, xerr=[_time_errm, _time_errp], yerr=_index_err,
                 fmt='o', label=sim_label)
    plt.axis([100., 1e6, None, None])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('grb130427_swift_pl_index', OUTPUT_FOLDER, False)
    """
    plt.figure('PL normalization')
    plt.errorbar(_phase, _norm, xerr=_phase_err, yerr=_norm_err, fmt='o',
                 label=sim_label)
    pl_normalization_spline.plot(show=False, label=mod_label)
    plt.axis([0., 1., None, None])
    plt.legend(bbox_to_anchor=(0.45, 0.95))
    if save:
        save_current_figure('crab_pl_norm', OUTPUT_FOLDER, False)
    """
    plt.show()


def run(save_plots=False):
    """Run all the tasks.
    """
    if os.path.exists(ANALYSIS_FILE_PATH):
        logger.info('%s exists, delete it if you want to recreate it.' %\
                    ANALYSIS_FILE_PATH)
    else:
        generate()
        prepare()
        analyze()
    plot(save_plots)


if __name__ == '__main__':
    run(save_plots=True)
