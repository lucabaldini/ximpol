#!/urs/bin/env python
import scipy as sp
from astropy import units as u

class xPsf():
    def __init__(self):
        self.psf=sp.random.normal
        pass
    def __call__(self,energy=0,theta=0):
        
        return self.psf(0.0,15.0)/3600.
