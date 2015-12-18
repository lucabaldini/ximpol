#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015                                                *
# * Nicola Omodei (nicola.omodei@stanford.edu)                        *
# * Melissa Pesce-Rollins (melissa.pesce.rollins@pi.infn.it)          *
# * Luca Baldini (luca.baldini@pi.infn.it)                            *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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



import time
import scipy
from scipy import interpolate
from matplotlib import pyplot as plt

from ximpol.srcmodel.xSource import xSource
from ximpol.srcmodel.xGenerator import xGenerator
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.irf.xAeff import xAeff
from ximpol.irf.xPsf import xPsf
from ximpol.irf.xModulation import xModulation
from ximpol.event.xEvent import xEvent
from ximpol.event.xEventList import xEventList
from ximpol.utils.xChrono import xChrono
from ximpol.__logging__ import logger, startmsg



def xobbsim(output_file_path, duration, start_time=0., time_steps=100,
            random_seed=0):
    """ ximpol fast simulator.
    """
    chrono = xChrono()
    logger.info('Loading the instrument response functions...')
    aeff       = xAeff()
    psf        = xPsf()
    modulation = xModulation()

    logger.info('Setting up the source model...')
    stop_time = start_time + duration
    emin=1
    emax=10
    phi0= 44.
    
    C=lambda t: 10.0*(1.0+scipy.cos(t))
    gamma=lambda t: -2.1
    
    mySource=xSource('Crab')
    ra0,dec0=mySource.getRADec()
    spectrum=xSpectralComponent('spectrum')
    
    times  = scipy.linspace(start_time, stop_time, time_steps)
    flux   = []
    events = []
    for t in times:
        spectrum.powerlaw(C(t),gamma(t))
        x,y = aeff.convolve(spectrum)
        f   = interpolate.UnivariateSpline(x,y,k=1,s=0)
        flux.append(f.integral(emin,emax))
        pass
    
    lc   = interpolate.UnivariateSpline(times,flux,k=1,s=0)
    logger.info('Done %s.' % chrono)
    logger.info('Extracting the event times...')
    S = xGenerator(lc,lc.integral)
    S.setMinMax(start_time,stop_time)    
    events_times = S.generate()
    logger.info('Done %s, %d events generated.' % (chrono, len(events_times)))

    event_list   = xEventList()
    logger.info('Entering the event loop...')
    for i,event_time in enumerate(events_times):
        _event      = xEvent()
        _event.time = event_time
        
        spectrum.powerlaw(C(event_time),gamma(event_time))
        energy_array,count_spectrum_array =  aeff.convolve(spectrum)
        
        f   = interpolate.UnivariateSpline(energy_array,count_spectrum_array,k=1,s=0)
        S   = xGenerator(f,f.integral)
        S.setMinMax(emin,emax)
        _event.energy=S.generate(1)[0]
        _ra,_dec=psf.smear(ra0,dec0)
        _event.setRADec(_ra,_dec)
        _event.angle= modulation.extract(_event.energy, phi0)
        event_list.fill(_event)        
        pass
    logger.info('Done %s.' % chrono)
    logger.info('Writing output file %s...' % output_file_path)
    event_list.write_fits(output_file_path)

    logger.info('Plotting stuff...')
    fig=plt.figure(figsize=(10,10),facecolor='w')
    #plt.subplots_adjust(hspace=0.001)
    ax = plt.subplot(221)
    plt.plot(times,flux)
    plt.xlabel('Time [s]')

    ax = plt.subplot(222)
    plt.hist2d(event_list.ra_array,event_list.dec_array,100,range=[[ra0-0.05,ra0+0.05],[dec0-0.05,dec0+0.05]])
    plt.xlabel('RA.')
    plt.ylabel('Dec.')

    ax = plt.subplot(223)
    plt.plot(event_list.time_array,event_list.energy_array,'o')
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [keV]')

    ax = plt.subplot(224)
    plt.hist(event_list.energy_array,bins=scipy.logspace(0,1,50),histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [keV]')
    logger.info('All done %s!' % chrono)
    plt.show()
    
    #spectrum.plot(scipy.linspace(1,10,100))
    #plt.xscale('log')
    #plt.yscale('log')
    #spectrum.polarization(0.1,89)
    #mySource.addComponent(spectrum)



if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        required=True,
                        help='the output FITS event file')
    parser.add_argument('-d', '--duration', type=float, default=None,
                        required=True,
                        help='the duration (in s) of the simulation')
    parser.add_argument('-t', '--start-time', type=float, default=0.,
                        help='the start time (MET in s) of the simulation')
    parser.add_argument('-n', '--time-steps', type=int, default=100,
                        help='the number of steps for sampling the lightcurve')
    parser.add_argument('-s', '--random-seed', type=int, default=0.,
                        help='the random seed for the simulation')
    args = parser.parse_args()
    startmsg()
    xobbsim(args.output_file, args.duration, args.start_time,
            args.time_steps, args.random_seed)
