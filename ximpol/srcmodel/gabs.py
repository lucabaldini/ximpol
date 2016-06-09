#!/urs/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
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
import numpy as np
import astropy.io.fits as pf
from astropy.coordinates import SkyCoord, Galactic
from astropy.wcs import WCS

from ximpol import XIMPOL_SRCMODEL
from ximpol.utils.logging_ import logger,abort

def get_column_density(ra, dec, HImap='LAB', lonlat = False):
    """Extract the column density from a HI map.

        Arguments
        ---------
        ra : float
            Right Ascension of the source [deg].

        dec: float
            Declination of the source [deg].

        map: str
            The label of the HI column density map to use. 
             Default: 'LAB' -> uses h1_nh_LAB.fits (LAB survay)
             Other options: 'DL' -> uses h1_nh_DL.fits (Dickey & Lockman)

        lonlat: Bool
             If True, ra and dec parameters are interpretated as galactic 
             L and B [deg]

        Warning
        -------
        Both the two HI column density maps available are expressed in galactic 
        coordinates. Hence, if RA and DEC are given, a convertion must be 
        performed.
    """
    nh_file_path = os.path.join(XIMPOL_SRCMODEL,'fits','h1_nh_%s.fits'%HImap)
    if not os.path.exists(nh_file_path):
        abort('ATT: HI column density maps %s not found!'%nh_file_path)
    if lonlat:
        logger.info('Getting HI column density for L = %f deg, B = %f deg, from %s HI map'\
                    %(ra,dec,HImap))
        column_density = nh_from_map(ra,dec,nh_file_path)
    else:
        logger.info('Getting HI column density for RA = %f deg, DEC = %f deg, from %s HI map'\
                    %(ra,dec,HImap))
        l,b = radec2lonlat(ra,dec)
        column_density = nh_from_map(l,b,nh_file_path)
    return column_density


def get_galactic_absorption(energy, column_density):
    """Returns the galactic transmission factor at a given energy and 
       column density.

        Arguments
        ---------
        energy : float
            Energy value [keV].

        column_density: float
             Column density at the source position
    """
    xsec = XsectionISM(energy)
    trans_factor = TransmissionFactor(xsec,column_density)
    return trans_factor
    

def radec2lonlat(ra, dec):
    """Function to get L and B if RA and DEC are given.
    """
    c = SkyCoord(ra,dec,unit='deg')
    cc = c.galactic
    cc.transform_to(Galactic)
    l,b = cc.l.degree, cc.b.degree
    return l,b


def nh_from_map(l, b, file_path):
    """Function to get NH value given RA and DEC and the HI map.
    """
    hdu = pf.open(file_path)
    w = WCS(hdu[0].header)
    row,col = w.wcs_world2pix(l,b,1)
    data = hdu[0].data
    column_density = data[int(col)][int(row)]
    return column_density


def XsectionISM(energy):
    """Returns the absorbtion x-section in cm^2 vs the energy in keV 
       Cross-section are defined between 0.03 keV and 10 keV
    """
    xsec_file = os.path.join(XIMPOL_SRCMODEL,'ascii','XsecFitParam.txt')
    if not os.path.exists(xsec_file):
        abort('ATT: cross-section fit parameters file %s not found!'%xsec_file)
    E_lo,E_hi,c0,c1,c2 = np.loadtxt(xsec_file,delimiter=',',unpack=True)
    bin_num = 0.
    for i in range(0,len(E_lo)):
        if E_lo[i] <= energy <= E_hi[i]:
            bin_num = i
            pass
        else:
            continue
    cross_section = (c0[bin_num]+c1[bin_num]*energy+c2[bin_num]*energy**2)\
               *energy**(-3)*10**(-24)
    return cross_section


def TransmissionFactor(cross_section, column_density):
    """Calculate the transmission factor for a given x-section and nh
    """
    factor = np.exp(-cross_section*column_density)
    return factor



def main():
    """
    """
    pass

if __name__ == '__main__':
    main()

