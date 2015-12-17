#!/usr/bin/env python
from astropy import coordinates
from astropy import units

class xEvent():
    def __init__(self):
        self.energy=0
        self.time=0
        self.ra=0
        self.dec=0
        self.x=0
        self.y=0
        self.l=0
        self.b=0
        self.angle=0
    def setRADec(self,ra,dec):
        self.ra=ra
        self.dec=dec
        #_coord=coordinates.SkyCoord(ra=ra*units.deg, dec=dec*units.deg, frame='fk5')
        #self.l=_coord.galactic.l.deg
        #self.b=_coord.galactic.b.deg

    def __str__(self):
        return 'Time=%10.3f Energy = %10.1f, Ra=%10.3f Dec=%10.3f Angle=%10.3f' %(self.time,self.energy,self.ra,self.dec,self.angle)
