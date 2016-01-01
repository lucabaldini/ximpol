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


import scipy
from scipy import integrate,optimize,interpolate
import sys

class xGenerator(object):
    def __init__(self,userFunction,userIntegral=None):
        self.userFunction         = userFunction
        if(userIntegral!=None):
            self.integral           = userIntegral
        else:
            self.integral           = lambda x,y: scipy.integrate.quad(self.userFunction,x,y,epsrel=1e-3)[0]
        pass

    def setMinMax(self,minX,maxX):
        self.minX          = minX
        self.maxX          = maxX
        self.totalIntegral = self.integral(self.minX,self.maxX)
        self.nEvents       = scipy.random.poisson(self.totalIntegral)

        def froot(u):
            f  = lambda x: self.integral(self.minX,x)/self.totalIntegral-u
            x0 = scipy.optimize.brentq(f,self.minX,self.maxX)
            return x0
        self.froot = froot
        pass

    def generate(self,nEvents=None):
        if nEvents is None: nEvents=self.nEvents
        X0                          = scipy.zeros(nEvents)
        j                           = 0
        uu                          = scipy.random.uniform(0,1,nEvents)
        for i in range(nEvents):
            X0[i]                     = self.froot(uu[i])
            pass
        #print("Done")
        X0.sort()
        self.events                 = X0
        return self.events

    def plot(self,nbins=300):
        from matplotlib import pyplot as plt
        counts,bin_edges          = scipy.histogram(self.events,nbins)
        bin_centres               = (bin_edges[:-1] + bin_edges[1:])/2.
        err                       = scipy.sqrt(counts)

        binsize                   = bin_edges[2]-bin_edges[1]
        xx                        = scipy.linspace(bin_edges[0],bin_edges[-1],1000)
        yy                        = scipy.array(self.userFunction(xx))*binsize
        #f,sub                     = plt.subplots(1,1,1)
        plt.plot(xx,yy,linewidth=1,color='#d95f02')
        plt.errorbar(bin_centres, counts, yerr=err, fmt='.',capsize=0,alpha=1,color='#7570b3')
        #plt.xlabel("Time (s)")
        #plt.ylabel("Counts")
        pass


def main():
    print 'all success!'


if __name__=='__main__':
    main()
