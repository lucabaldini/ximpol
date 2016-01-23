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
from ximpol.utils.logging_ import logger, startmsg
from ximpol import XIMPOL_IRF
from ximpol.srcmodel.spectrum import xCountSpectrum


def xpobssim(output_file_path, config_file_path, irf_name, duration, start_time,
             time_steps, random_seed):
    """Run the ximpol fast simulator.
    """
    assert(config_file_path.endswith('.py'))
    if output_file_path is None:
        file_name = os.path.basename(config_file_path).replace('.py', '.fits')
        output_file_path = file_name

    chrono = xChrono()
    logger.info('Setting the random seed to %d...' % random_seed)
    numpy.random.seed(random_seed)

    logger.info('Loading the instrument response functions...')
    aeff, psf, modf, edisp = load_irfs(irf_name)
    logger.info('Done %s.' % chrono)

    logger.info('Setting up the source model...')
    module_name = os.path.basename(config_file_path).replace('.py', '')
    ROI_MODEL = imp.load_source(module_name, config_file_path).ROI_MODEL
    stop_time = start_time + duration
    sampling_time = numpy.linspace(start_time, stop_time, time_steps)
    event_list=ROI_MODEL.rvs_event_list( aeff, psf, modf, edisp,sampling_time)
    logger.info('Done %s.' % chrono)
    event_list.write_fits(output_file_path)
    logger.info('All done %s!' % chrono)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        help='the output FITS event file')
    parser.add_argument('-c', '--config-file', type=str, required=True,
                        help='the input configuration file')
    parser.add_argument('-r', '--irf-name', type=str, default='xipe_baseline',
                        help='the input configuration file')
    parser.add_argument('-d', '--duration', type=float, default=10,
                        help='the duration (in s) of the simulation')
    parser.add_argument('-t', '--start-time', type=float, default=0.,
                        help='the start time (MET in s) of the simulation')
    parser.add_argument('-n', '--time-steps', type=int, default=100,
                        help='the number of steps for sampling the lightcurve')
    parser.add_argument('-s', '--random-seed', type=int, default=0,
                        help='the random seed for the simulation')
    args = parser.parse_args()
    startmsg()
    xpobssim(args.output_file, args.config_file, args.irf_name, args.duration,
             args.start_time, args.time_steps, args.random_seed)
