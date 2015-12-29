#!/urs/bin/env python
import scipy as sp
from astropy import units as u

class xPsf():
    def __init__(self,rms_arcsec=15.):
        self.psf = sp.random.normal
        self.rms = rms_arcsec/3600.
        
        pass
    def __call__(self,energy=0,theta=0):
        
        return self.psf(0.0,self.rms)
    
    def smear(self,ra0,dec0):
        radius=self()
        angle=sp.random.uniform(0,2*sp.pi)
        return  ra0+radius*sp.cos(angle),dec0+radius*sp.sin(angle),
