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
from matplotlib import pyplot as plt

from ximpol.srcmodel.xSource import xSource
from ximpol.srcmodel.xGenerator import xGenerator
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.irf.xAeff import xAeff
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


def xpobssim(output_file_path, duration, start_time, time_steps, random_seed):
    """Run the ximpol fast simulator.
    """
    chrono = xChrono()
    logger.info('Setting the random seed to %d...' % random_seed)
    numpy.random.seed(random_seed)
    logger.info('Loading the instrument response functions...')
    aeff = xAeff()
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.psf')
    psf = xPointSpreadFunction(file_path)
    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    modf = xModulationFactor(file_path)
    logger.info('Done %s.' % chrono)
    logger.info('Setting up the source model...')

    stop_time = start_time + duration
    emin=1
    emax=10
    phi0= 44.

    # This is still needed for the spectrum.powerlaw(C(t),gamma(t))
    # but eventually should be removed.
    C = lambda t: 10.0*(1.0 + numpy.cos(t))
    gamma = lambda t: -2.0 + 0.01*t

    mySource = xSource('Crab', resolve_name=False)
    ra0, dec0 = mySource.getRADec()
    spectrum = xSpectralComponent('spectrum')


    # Set up the generator for the spectrum.
    # Note we accidentally left out the convolution with the effective
    # area, need to put that back in!!!!

    def dNdE(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    E = numpy.linspace(1, 10, 100)
    t = numpy.linspace(0, 100, 100)
    fmt = dict(auxname='Time', auxunits='s', rvname='Energy', rvunits='keV',
               pdfname='dN/dE')
    spec_gen = xUnivariateAuxGenerator(t, E, dNdE, **fmt)

    times  = numpy.linspace(start_time, stop_time, time_steps)
    flux   = []
    events = []
    for t in times:
        spectrum.powerlaw(C(t),gamma(t))
        x, y = aeff.convolve(spectrum)
        f = interpolate.UnivariateSpline(x,y,k=1,s=0)
        flux.append(f.integral(emin,emax))
        pass

    lc = interpolate.UnivariateSpline(times,flux,k=1,s=0)
    logger.info('Done %s.' % chrono)
    logger.info('Extracting the event times...')
    S = xGenerator(lc,lc.integral)
    S.setMinMax(start_time,stop_time)
    events_times = S.generate()
    logger.info('Done %s, %d events generated.' % (chrono, len(events_times)))

    event_list = xMonteCarloEventList()
    event_list.set_column('TIME', events_times)
    event_list.set_column('ENERGY', spec_gen.rvs(events_times))
    ra, dec = psf.smear_single(ra0, dec0, len(event_list))
    event_list.set_column('RA', ra)
    event_list.set_column('DEC', dec)

    # And this seems to be fundamentally different from the
    # xUnivariateAuxGenerator case, so an intermediate layer might be
    # necessary.
    _x = modf.x.copy()
    _y = numpy.linspace(0, 2*numpy.pi, 100)
    _z = numpy.zeros(shape = (_x.size, _y.size))
    for i, _xp in enumerate(_x):
        mu = modf(_xp)
        for j, _yp in enumerate(_y):
            _z[i, j] = (1.0 - mu)/2*numpy.pi +\
                       mu/numpy.pi*numpy.power(numpy.cos(_yp - phi0), 2.)

    modf_spline = xInterpolatedBivariateSplineLinear(_x, _y, _z)
    modf_ppf = modf_spline.build_vppf()

    pe_angles = modf_ppf(event_list['ENERGY'],\
                         numpy.random.sample(len(event_list)))
    event_list.set_column('PE_ANGLE', pe_angles)
    logger.info('Done %s.' % chrono)
    logger.info('Writing output file %s...' % output_file_path)
    event_list.write_fits(output_file_path)

    logger.info('Plotting stuff...')
    fig=plt.figure(figsize=(10,10), facecolor='w')
    ax = plt.subplot(221)
    plt.plot(times, flux)
    plt.xlabel('Time [s]')

    ax = plt.subplot(222)
    fov = 50./3600
    plt.hist2d(event_list['RA'], event_list['DEC'],100,
               range=[[ra0-fov,ra0+fov],[dec0-fov,dec0+fov]])
    plt.xlabel('RA.')
    plt.ylabel('Dec.')

    ax = plt.subplot(223)
    plt.plot(event_list['TIME'], event_list['ENERGY'], 'o')
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [keV]')

    ax = plt.subplot(224)
    plt.hist(event_list['ENERGY'], bins=numpy.logspace(0,1,50),histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [keV]')
    logger.info('All done %s!' % chrono)
    plt.show()


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
