#!/usr/bin/env python
from astropy import units as u
from astropy import coordinates as coord
import scipy as sp

class xSource():
    '''Source class:
    name of the source can be automatically resolved by NED and SIMBAD
    astropy.coordinates handle coordinate transformations
    astropy.units define units for the different quantities
    
    author: nicola.omodei@stanford.edu
    '''
    def __init__(self,name,resolve_name=True):
        self.name=name
        if resolve_name: 
            try:
                self.coord=coord.SkyCoord.from_name(name)
            except:
                print 'WARNING [xSource]: Unable to find coordinates for name %s' % name 
                self.setLB(0,0)  
            pass
        else: self.setLB(0,0)
        self.components=sp.array([])
        
    def setSkyCord(self,coord):
        self.coord=coord

    def setRADec(self,ra,dec):
        self.setSkyCord(coord.SkyCoord(ra=ra*u.deg, dec=dec*u.deg, frame='fk5'))
    
    def setLB(self,l,b):
        self.setSkyCord(coord.SkyCoord(l=l*u.deg, b=b*u.deg, frame='galactic'))

    def getSetSkyCord(self): 
        return self.coord
    
    def getRADec(self): 
        return self.coord.fk5.ra.deg,self.coord.fk5.dec.deg
    
    def getLB(self): 
        return self.coord.galactic.l.deg,self.coord.galactic.b.deg
    
    def addComponent(self,component):
        print 'Component %s added' % component.name
        self.components=sp.append(self.components,component)
        pass

def test():
    print 'All Done!'

if __name__=='__main__':
    test()
