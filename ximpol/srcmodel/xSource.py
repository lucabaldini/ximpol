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


from astropy import units as u
from astropy import coordinates as coord
import scipy as sp


class xSource():
    """Source class:
    name of the source can be automatically resolved by NED and SIMBAD
    astropy.coordinates handle coordinate transformations
    astropy.units define units for the different quantities

    author: nicola.omodei@stanford.edu
    """
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


def main():
    print 'All Done!'


if __name__=='__main__':
    main()
