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
from scipy import interpolate

from ximpol.srcmodel.xSource import xSource
from ximpol.srcmodel.xGenerator import xGenerator
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.irf.xAeff import xAeff
from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.psf import xPointSpreadFunction
from ximpol.irf.mrf import xModulationFactor
from ximpol.event import xMonteCarloEventList
from ximpol.utils.profile import xChrono
from ximpol.utils.logging_ import logger, startmsg
from ximpol import XIMPOL_IRF
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.core.rand import xUnivariateAuxGenerator
from ximpol.srcmodel.spectrum import xCountSpectrum


def xpobssim(output_file_path, duration, start_time, time_steps, random_seed):
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
    logger.info('Done %s.' % chrono)

    logger.info('Setting up the source model...')
    stop_time = start_time + duration
    min_energy = 1
    max_energy = 10
    polarization_angle = 44.
    polarization_degree = 1.

    mySource = xSource('Crab', resolve_name=False)
    source_ra, source_dec = mySource.getRADec()

    def dNdE(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    t = numpy.linspace(start_time, stop_time, time_steps)
    count_spectrum = xCountSpectrum(dNdE, aeff, t)
    modf.build_generator(polarization_angle, polarization_degree)
    logger.info('Done %s.' % chrono)

    logger.info('Extracting the event times...')
    num_events = numpy.random.poisson(count_spectrum.light_curve.norm())
    event_times = count_spectrum.light_curve.rvs(num_events)
    logger.info('Done %s, %d events generated.' % (chrono, num_events))

    logger.info('Filling output columns...')
    event_list = xMonteCarloEventList()
    event_list.set_column('TIME', event_times)
    event_list.set_column('ENERGY', count_spectrum.rvs(event_times))
    ra, dec = psf.smear_single(source_ra, source_dec, num_events)
    event_list.set_column('RA', ra)
    event_list.set_column('DEC', dec)
    event_list.set_column('PE_ANGLE', modf.rvs(event_list['ENERGY']))
    logger.info('Done %s.' % chrono)

    event_list.write_fits(output_file_path)
    logger.info('All done %s!' % chrono)


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        required=True,
                        help='the output FITS event file')
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
    xpobssim(args.output_file, args.duration, args.start_time,
            args.time_steps, args.random_seed)
