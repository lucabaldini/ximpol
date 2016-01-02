#!/usr/bin/env python
import scipy as sp
from matplotlib import pyplot as plt
from ximpol.srcmodel.xSpectralComponent import xSpectralComponent
from ximpol.srcmodel.xSource import xSource

# EXAMPLE:
Crab=xSource('Crab')#,resolve_name=False)
#Crab.setRADec(83.63309062468973, 22.014494786634714)
pl1 = xSpectralComponent('pl1')
pl2 = xSpectralComponent('pl2')
co1 = xSpectralComponent('co1')

pl1.powerlaw(1.0,-2)
pl2.powerlaw(1.0,-1)

pl3=pl1+pl2

co1.highecut(2.0,5.0)

synch=(pl1+pl2)*co1
synch.polarization(degree=0.17,angle=67)

Crab.addComponent(synch)

xx=sp.linspace(1,10,100)
pl1.plot(xx)
pl2.plot(xx)
pl3.plot(xx,':')
synch.plot(xx,'r')

plt.xscale('log')
plt.yscale('log')
plt.legend()
plt.show()

#Crab.AddComponent(...)
#Crab.TemporalProfile()

#CygX1=xSource('Cyg X-1')
#Crab.SetRADec(83.633083,22.014500)
#print Crab.getLB()
#print Crab.getRADec()
#l,b=Crab.getLB()
#Crab.setLB(l,b)
#print Crab.getLB()
#print Crab.getRADec()
#Angle=Crab.coord.separation(CygX1.coord)
#print Angle.deg
