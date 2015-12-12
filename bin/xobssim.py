#!/usr/bin/env python
from ximpol.srcmodel.xSource import xSource
from ximpol.srcmodel.xGenerator import xGenerator
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.irf.xAeff import xAeff
from ximpol.irf.xPsf import xPsf

import scipy as sp  
from scipy import interpolate
from matplotlib import pyplot as plt

if __name__=='__main__':
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
    
    times=sp.linspace(tmin,tmax,100)
    flux=[]
    events=[]
    for t in times:
        spectrum.powerlaw(C(t),gamma(t))
        x,y = aeff.convolve(spectrum)
        f   = interpolate.UnivariateSpline(x,y,k=1,s=0)
        flux.append(f.integral(emin,emax))
        pass
    
    lc   = interpolate.UnivariateSpline(times,flux,k=1,s=0)
    
    S=xGenerator(lc,lc.integral)
    S.setMinMax(tmin,tmax)
    events_times=S.generate()
    events_energies=[]
    events_position_ra=[]
    events_position_dec=[]
    events_photoelectron_angle=[]
    
    print 'generating %d events...' % len(events_times)
    for i,evt in enumerate(events_times):
        #print '%.2f' % (i/(1.0*len(events_times)))
        spectrum.powerlaw(C(evt),gamma(evt))
        x,y = aeff.convolve(spectrum)
        f   = interpolate.UnivariateSpline(x,y,k=1,s=0)
        S   = xGenerator(f,f.integral)
        S.setMinMax(emin,emax)
        events_energies.append(S.generate(1)[0])
        #ra,dec=mySource.getRADec()
        
        events_position_ra.append(ra0+psf())
        events_position_dec.append(dec0+psf())
        events_photoelectron_angle.append(0.0)        
        pass
    
    fig=plt.figure(figsize=(10,10),facecolor='w')
    #plt.subplots_adjust(hspace=0.001)
    ax = plt.subplot(221)
    plt.plot(times,flux)
    plt.xlabel('Time [s]')

    ax = plt.subplot(222)
    plt.hist2d(events_position_ra,events_position_dec,100,range=[[ra0-0.05,ra0+0.05],[dec0-0.05,dec0+0.05]])
    plt.xlabel('RA.')
    plt.ylabel('Dec.')

    ax = plt.subplot(223)
    plt.plot(events_times,events_energies,'o')
    plt.yscale('log')
    plt.xlabel('Time [s]')
    plt.ylabel('Energy [keV]')

    ax = plt.subplot(224)
    plt.hist(events_energies,bins=sp.logspace(0,1,50),histtype='step')
    plt.xscale('log')
    plt.yscale('log')
    plt.xlabel('Energy [keV]')

    plt.show()
    
    #spectrum.plot(sp.linspace(1,10,100))
    #plt.xscale('log')
    #plt.yscale('log')
    #spectrum.polarization(0.1,89)
    #mySource.addComponent(spectrum)
