#!/urs/bin/env python
import scipy as sp
from astropy import units as u
from astropy import coordinates as coord

class xSpectralComponent():
    '''Spectral Component class:
    define the spectral component that can be used to build a model spectrum to use with the simulator
    author: nicola.omodei@stanford.edu
    '''
    def __init__(self,name):
        self.name=name
        self.degree=0
        self.angle=0
        self.f=None
        pass
    
    def powerlaw(self,C,gamma):
        'C is the normalization at 1keV'
        self.f=lambda x: C*sp.power(x,gamma)
        pass

    def highecut(self,foldE,cutoffE=0.0):
        'C is the normalization at 1keV'
        self.f=lambda x: sp.exp(-sp.maximum(0.0,(x-cutoffE)/foldE))
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

def test():
    print 'all success!'

if __name__=='__main__':
    test()
