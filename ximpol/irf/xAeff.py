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


class xAeff():
    def __init__(self,myFile=os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')):
        hdulist = fits.open(myFile)
        data=hdulist['SPECRESP'].data
        self.ENERG_LO=data.field('ENERG_LO')
        self.ENERG_HI=data.field('ENERG_HI')
        self.SPECRESP=data.field('SPECRESP') # cm^2
        self.CentralEnergy=(self.ENERG_HI+self.ENERG_LO)*0.5
        pass

    def convolve(self,dNdE):
        # dNdE is in ph/cm^2/s/keV
        return self.CentralEnergy,self.SPECRESP*dNdE(self.CentralEnergy)

    def __call__(self,x):
        pass

    def plot(self,**kwargs):
        from matplotlib import pyplot as plt
        plt.plot(self.CentralEnergy,self.SPECRESP,**kwargs)
        plt.show()


def main():
    aeff=xAeff()
    aeff.plot()

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
