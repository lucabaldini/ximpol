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


import os
import astropy
from astropy.io import fits
from ximpol import XIMPOL_IRF
import scipy as sp
from scipy import interpolate
from ximpol.srcmodel.xGenerator import xGenerator


class xModulation():
    def __init__(self,myFile=os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')):
        hdulist = fits.open(myFile)
        data=hdulist['MODFRESP'].data
        self.ENERG_LO=data.field('ENERG_LO')
        self.ENERG_HI=data.field('ENERG_HI')
        self.MODFRESP=data.field('MODFRESP')
        self.CentralEnergy=(self.ENERG_HI+self.ENERG_LO)*0.5
        self.modulation_spline = interpolate.UnivariateSpline(self.CentralEnergy,self.MODFRESP,k=1,s=0)

        pass

    def __call__(self,x):
        pass

    def plot(self,**kwargs):
        from matplotlib import pyplot as plt
        plt.plot(self.CentralEnergy,self.MODFRESP,**kwargs)
        plt.show()

    def extract(self,energy, phi0):
        mu  = self.modulation_spline(energy)
        phis=sp.linspace(0,2*sp.pi,100)
        muf = (1.0-mu)/2*sp.pi+mu/sp.pi*sp.power(sp.cos(phis-phi0),2.)
        f=interpolate.UnivariateSpline(phis,muf,k=1,s=0)
        S   = xGenerator(f,f.integral)
        S.setMinMax(0,2*sp.pi)
        return S.generate(1)[0]


def main():
    from scipy import random
    modulation=xModulation()
    modulation.plot()
    nevents=1000
    for i in range(nevents):
        energy=random.uniform(1,10)
        print modulation.extract(energy,76)
    #x,y=aeff.convolve(Crab.components[0])
    #f = interpolate.UnivariateSpline(x,y,k=1,s=0)
    #plt.plot(x,Crab.components[0](x))

    #plt.plot(x,y)
    #plt.plot(x,f(x))
    #S=xSimulator(f,f.integral)
    #S.setMinMax(1,10)

    #S.generate()
    ##S.plot()
    #plt.xscale('log')
    #plt.yscale('log')
    #S.nEvents

if __name__=='__main__':
    main()
