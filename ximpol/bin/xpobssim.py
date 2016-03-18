#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
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


__description__ = 'Run the ximpol fast simulator'


import os
import numpy
import imp

from ximpol.irf import load_irfs
from ximpol.evt.event import xMonteCarloEventList
from ximpol.utils.profile import xChrono
from ximpol.utils.os_ import mkdir
from ximpol.utils.logging_ import logger, startmsg
from ximpol import XIMPOL_IRF, XIMPOL_DATA
from ximpol.srcmodel.spectrum import xCountSpectrum


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__,
                                 formatter_class=formatter)
PARSER.add_argument('--outfile', type=str, default=None,
                    help='the output FITS event file')
PARSER.add_argument('--configfile', type=str, required=True,
                    help='the input configuration file')
PARSER.add_argument('--irfname', type=str, default='xipe_baseline',
                    help='the input configuration file')
PARSER.add_argument('--duration', type=float, default=10,
                    help='the duration (in s) of the simulation')
PARSER.add_argument('--tstart', type=float, default=0.,
                    help='the start time (MET in s) of the simulation')
PARSER.add_argument('--seed', type=int, default=0,
                    help='the random seed for the simulation')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')


class xSimulationInfo:

    """Empty container to pass along all the relevant information about the
    simulation.
    """

    pass


def xpobssim(**kwargs):
    """Run the ximpol fast simulator.
    """
    assert(kwargs['configfile'].endswith('.py'))
    if kwargs['outfile'] is None:
        outfile = os.path.basename(kwargs['configfile']).replace('.py', '.fits')
        mkdir(XIMPOL_DATA)
        kwargs['outfile'] = os.path.join(XIMPOL_DATA, outfile)
        logger.info('Setting output file path to %s...' % kwargs['outfile'])
    if os.path.exists(kwargs['outfile']) and not kwargs['clobber']:
        logger.info('Output file %s already exists.' % kwargs['outfile'])
        logger.info('Remove the file or set "clobber = True" to overwite it.')
        return kwargs['outfile']
    chrono = xChrono()
    logger.info('Setting the random seed to %d...' % kwargs['seed'])
    numpy.random.seed(kwargs['seed'])
    logger.info('Loading the instrument response functions...')
    aeff, psf, modf, edisp = load_irfs(kwargs['irfname'])
    logger.info('Done %s.' % chrono)
    logger.info('Setting up the source model...')
    module_name = os.path.basename(kwargs['configfile']).replace('.py', '')
    ROI_MODEL = imp.load_source(module_name, kwargs['configfile']).ROI_MODEL
    logger.info(ROI_MODEL)
    if kwargs['tstart'] < ROI_MODEL.min_validity_time():
        kwargs['tstart'] = ROI_MODEL.min_validity_time()
        logger.info('Simulation start time set to %s...' % kwargs['tstart'])
    tstop = kwargs['tstart'] + kwargs['duration']
    if tstop > ROI_MODEL.max_validity_time():
        tstop = ROI_MODEL.max_validity_time()
        logger.info('Simulation stop time set to %s...' % tstop)
    gti_list = [(kwargs['tstart'], tstop)]
    kwargs['tstop'] = tstop
    event_list = ROI_MODEL.rvs_event_list(aeff, psf, modf, edisp, **kwargs)
    logger.info('Done %s.' % chrono)
    simulation_info = xSimulationInfo()
    simulation_info.gti_list = gti_list
    simulation_info.roi_model = ROI_MODEL
    simulation_info.irf_name = kwargs['irfname']
    simulation_info.aeff = aeff
    simulation_info.psf = psf
    simulation_info.modf = modf
    simulation_info.edisp = edisp
    event_list.write_fits(kwargs['outfile'], simulation_info)
    logger.info('All done %s!' % chrono)
    return kwargs['outfile']


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpobssim(**args.__dict__)
