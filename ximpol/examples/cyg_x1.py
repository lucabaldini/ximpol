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
from ximpol.config.cyg_x1 import PL_NORM, PL_INDEX, pol_degree_spline, pol_angle_spline

"""Script-wide simulation and analysis settings.
"""
CFG_FILE = os.path.join(XIMPOL_CONFIG, 'cyg_x1.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'cyg_x1')
MCUBE_FILE_PATH = os.path.join(XIMPOL_DATA, 'cyg_x1_mcube.fits')
ANALYSIS_FILE_PATH = '%s_analysis.txt' % OUT_FILE_PATH_BASE
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
SIM_DURATION = 50000.


E_BINNING = [2.0, 4.0, 6.0, 8.0]
"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)


def run():
    #First simulate the events
    PIPELINE.xpobssim(configfile=CFG_FILE, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)
    #Bin into PHA1 format for xspec
    #pha1_file_path = PIPELINE.xpbin(EVT_FILE_PATH, algorithm='PHA1')
    #spec_fitter = PIPELINE.xpxspec(pha1_file_path,emin=1.0)
    #Fit with xspec
    #(index, index_err), (norm, norm_err) = spec_fitter.fit_parameters()
    #logger.info('Fitted PL norm = %.4f +- %4f (input = %.4f)' %\
    #            (norm, norm_err,PL_NORM))
    #logger.info('Fitted PL index = %.4f +- %4f (input = %.4f)' %\
    #            (index, index_err, PL_INDEX))
    #make the mcube
    PIPELINE.xpbin(EVT_FILE_PATH, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING)

def analyze():

    """Analyze the data.Testing this method, but I must be missing something, it does not work yet.
    """
    logger.info('Opening output file %s...' % ANALYSIS_FILE_PATH)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
    _mcube = xBinnedModulationCube(MCUBE_FILE_PATH)
    _mcube.fit()
    _fit_results = _mcube.fit_results[0]
    _pol_deg = _fit_results.polarization_degree
    _pol_deg_err = _fit_results.polarization_degree_error
    _pol_angle = _fit_results.phase
    _pol_angle_err = _fit_results.phase_error
    #_energy_mean = _mcube.emean[0]
    #_emin = _mcube.emean[1]
    #_emax = _mcube.emean[2]
    #_data = (_energy_mean,_emin, _emax,_pol_deg, _pol_deg_err, _pol_angle, _pol_angle_err)
    print _data
    _fmt = ('%.4e   ' * len(_data)).strip()
    _fmt = '%s\n' % _fmt
    _line = _fmt % _data
    analysis_file.write(_line)
    analysis_file.close()



def view():
    #_energy_mean,_emin, _emax, _pol_deg, _pol_deg_err, _pol_angle, \
    #    _pol_angle_err = \
    #                    numpy.loadtxt(ANALYSIS_FILE_PATH, unpack=True)
    _mcube = xBinnedModulationCube(MCUBE_FILE_PATH)
    _mcube.fit()
    _fit_results = _mcube.fit_results[0]
    plt.figure('Polarization degree')
    _mcube.plot_polarization_degree(show=False, color='blue')
    pol_degree_spline.plot(color='lightgray',label='Model', show=False)

    #plt.errorbar(_energy_mean, _pol_deg, yerr=_pol_deg_err, color='blue',marker='o')


    plt.figure('Polarization angle')
    _mcube.plot_polarization_angle(show=False, color='blue', degree=False)
    pol_angle_spline.plot(color='lightgray',label='Model', show=False)
    #plt.errorbar(_energy_mean,_pol_angle, yerr= _pol_angle_err,color='blue',marker='o')

    plt.legend()
    plt.show()

if __name__ == '__main__':
    run()
    #analyze()
    view()
