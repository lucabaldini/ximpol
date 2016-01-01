#!/urs/bin/env python
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


import scipy


class xSpectralComponent():
    """Spectral Component class:
    define the spectral component that can be used to build a model spectrum to use with the simulator
    author: nicola.omodei@stanford.edu
    """
    def __init__(self,name):
        self.name=name
        self.degree=0
        self.angle=0
        self.f=None
        pass

    def powerlaw(self,C,gamma):
        'C is the normalization at 1keV'
        self.f=lambda x: C*scipy.power(x,gamma)
        pass

    def highecut(self,foldE,cutoffE=0.0):
        'C is the normalization at 1keV'
        self.f=lambda x: scipy.exp(-scipy.maximum(0.0,(x-cutoffE)/foldE))
        pass



    def __call__(self,x):
        return self.f(x)

    def __add__(self,other):

        newname='%s+%s'%(self.name,other.name)
        summed=xSpectralComponent(newname)
        summed.f=lambda x: self.f(x)+other.f(x)
        return summed

    def __mul__(self,other):
        newname='(%s)*(%s)'%(self.name,other.name)
        multiplied=xSpectralComponent(newname)
        multiplied.f=lambda x: self.f(x)*other.f(x)
        return multiplied

    def polarization(self,degree,angle):
        self.degree=degree
        self.angle =angle

    def plot(self,x,*args, **kwargs):
        from matplotlib import pyplot as plt
        plt.plot(x,self.f(x),label=self.name,*args, **kwargs)
        pass


def main():
    print 'all success!'


if __name__=='__main__':
    main()
