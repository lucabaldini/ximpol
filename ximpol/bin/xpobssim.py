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

from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.psf import xPointSpreadFunction
from ximpol.irf.mrf import xModulationFactor
from ximpol.irf.rmf import xEnergyDispersion
from ximpol.evt.event import xMonteCarloEventList
from ximpol.utils.profile import xChrono
from ximpol.utils.logging_ import logger, startmsg
from ximpol import XIMPOL_IRF
from ximpol.srcmodel.spectrum import xCountSpectrum


def xpobssim(output_file_path, config_file_path, duration, start_time,
             time_steps, random_seed):
    """Run the ximpol fast simulator.
    """
    chrono = xChrono()
    logger.info('Setting the random seed to %d...' % random_seed)
    numpy.random.seed(random_seed)

    logger.info('Loading the instrument response functions...')
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
    aeff = xEffectiveArea(file_path)
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    psf = xPointSpreadFunction(file_path)
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    modf = xModulationFactor(file_path)
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.rmf')
    edisp = xEnergyDispersion(file_path)
    logger.info('Done %s.' % chrono)

    logger.info('Setting up the source model...')
    module_name = os.path.basename(config_file_path).replace('.py', '')
    source = imp.load_source(module_name, config_file_path).source
    stop_time = start_time + duration
    t = numpy.linspace(start_time, stop_time, time_steps)
    count_spectrum = xCountSpectrum(source.spectrum, aeff, t)

    modf.build_generator(source.polarization_angle, source.polarization_degree)
    logger.info('Done %s.' % chrono)

    logger.info('Extracting the event times...')
    event_list = xMonteCarloEventList()
    num_events = numpy.random.poisson(count_spectrum.light_curve.norm())
    _time = count_spectrum.light_curve.rvs(num_events)
    logger.info('Sorting event times...')
    _time.sort()
    logger.info('Done %s, %d events generated.' % (chrono, num_events))
    logger.info('Filling output columns...')
    event_list.set_column('TIME', _time)
    _mc_energy = count_spectrum.rvs(_time)
    event_list.set_column('MC_ENERGY', _mc_energy)
    _pha = edisp.matrix.rvs(_mc_energy)
    event_list.set_column('PHA', _pha)
    event_list.set_column('ENERGY', edisp.ebounds(_pha))
    _ra, _dec = psf.smear_single(source.ra, source.dec, num_events)
    event_list.set_column('RA', _ra)
    event_list.set_column('DEC', _dec)
    _pe_angle = modf.rvs(_mc_energy)
    event_list.set_column('PE_ANGLE', _pe_angle)
    event_list.set_column('MC_RA', source.ra)
    event_list.set_column('MC_RA', source.dec)
    event_list.set_column('MC_SRC_ID', source.identifier)
    logger.info('Done %s.' % chrono)
    event_list.write_fits(output_file_path)
    logger.info('All done %s!' % chrono)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-o', '--output-file', type=str, required=True,
                        help='the output FITS event file')
    parser.add_argument('-c', '--config-file', type=str, required=True,
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
    xpobssim(args.output_file, args.config_file, args.duration, args.start_time,
            args.time_steps, args.random_seed)
