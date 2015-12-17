#!/usr/bin/env python
import time

from ximpol.srcmodel.xSource import xSource
from ximpol.srcmodel.xGenerator import xGenerator
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.irf.xAeff import xAeff
from ximpol.irf.xPsf import xPsf
from ximpol.event.xEvent import xEvent
from ximpol.event.xEventList import xEventList

import scipy as sp  
from scipy import interpolate
from matplotlib import pyplot as plt

if __name__=='__main__':
    t0=time.time()
    tmin=0
    tmax=10
    emin=1
    emax=10
    
    aeff=xAeff()
    psf=xPsf()

    C=lambda t: 10.0*(1.0+sp.cos(t))
    gamma=lambda t: -2.1
    
    mySource=xSource('Crab')
    ra0,dec0=mySource.getRADec()

    
    print mySource.getLB()
    print ra0, dec0
    spectrum=xSpectralComponent('spectrum')
    
    times  = sp.linspace(tmin,tmax,100)
    flux   = []
    events = []
    for t in times:
        spectrum.powerlaw(C(t),gamma(t))
        x,y = aeff.convolve(spectrum)
        f   = interpolate.UnivariateSpline(x,y,k=1,s=0)
        flux.append(f.integral(emin,emax))
        pass
    
    lc   = interpolate.UnivariateSpline(times,flux,k=1,s=0)
    
    S = xGenerator(lc,lc.integral)
    S.setMinMax(tmin,tmax)
    
    events_times = S.generate()
    event_list   = xEventList()
    
    print 'generating %d events...' % len(events_times)
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
        _event.angle=0.0
        event_list.fill(_event)        
        pass
    event_list.write_fits('test.fits')
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
    plt.hist(event_list.energy_array,bins=sp.logspace(0,1,50),histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [keV]')
    print time.time()-t0
    plt.show()
    
    #spectrum.plot(sp.linspace(1,10,100))
    #plt.xscale('log')
    #plt.yscale('log')
    #spectrum.polarization(0.1,89)
    #mySource.addComponent(spectrum)
