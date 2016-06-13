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
import numpy
import astropy.io.fits as pf
from astropy.coordinates import SkyCoord, Galactic
from astropy.wcs import WCS

from ximpol import XIMPOL_SRCMODEL
from ximpol.utils.logging_ import logger,abort
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


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
    E_lo,E_hi,c0,c1,c2 = numpy.loadtxt(xsec_file,delimiter=',',unpack=True)
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
    factor = numpy.exp(-cross_section*column_density)
    return factor



class xpeInterstellarAbsorptionModel:

    """Class implementing the insterstellar absorption model using the
    Wisconsin (Morrison and McCammon; ApJ 270, 119) cross-sections.

    Here we essentially read the numbers in table 2 from the paper and
    build an interpolated univariate linear spline with the photoabsorption
    cross section values.
    """

    def __init__(self, num_samples=250):
        """Constructor.

        Arguments
        ---------
        num_samples : int
            The number of data points used to sample (logarithmically) the
            photon energies when tabulating the absorption cross section
            internally.
        """
        file_path = os.path.join(XIMPOL_SRCMODEL, 'ascii', 'XsecFitParam.txt')
        if not os.path.exists(file_path):
            abort('Could not find %s' % file_path)
        # Read the data from file and calculate the minimum and maximum
        # energies.
        e_lo, e_hi, c0, c1, c2 = numpy.loadtxt(file_path, delimiter=',',
                                               unpack=True)
        emin = e_lo[0]
        emax = e_hi[-1]
        # Sample the energy logarithmically between emin and emax.
        _x = numpy.logspace(numpy.log10(emin), numpy.log10(emax), num_samples)
        # Here comes the interesting part---the cross section is tabulated
        # by means of a set of piecewise quadratic functions, and the
        # three coefficients are provided in each energy bin. We do some
        # magic with the numpy.digitize() function, returning for each value
        # in _x its bin index in the original table. Note that we need to
        # clip the array removing negative numbers, as the double log/exp
        # operations done in defining the binning where driving the first
        # index to -1.
        _bin = (numpy.digitize(_x, e_lo) - 1).clip(min=0)
        # Calculate the cross section values---compact, isn't it?
        _y = 1.0e-24*(c0[_bin] + c1[_bin]*_x + c2[_bin]*(_x**2.))/(_x**3.)
        # And, finally, build the basic spline underlying the model.
        _fmt = dict(xname='Energy', xunits='keV',
                    yname='ISM absorption cross section', yunits='cm$^2$')
        self.xsection = xInterpolatedUnivariateSplineLinear(_x, _y, **_fmt)

    def __call__(self, energy):
        """Return the value(s) of the photoelectric absorption cross section
        of the model for a value (or an array of values) of energy.

        Arguments
        ---------
        energy : float or array
            The energy (or energies) at which the cross section should be
            calculated.
        """
        return self.xsection(energy)

    def transmission_factor(self, column_density):
        """Return the transmission factor for a given column density.

        This is essentially returning

        .. math::
            \\varepsilon = \\exp(-n_H\\sigma)

        Arguments
        ---------
        column_density : float
            The column density at which the transmission factor should be
            calculated.

        Warning
        -------
        We do have an issue with the extrapolation, here, as at the current
        stage there is no guarantee that the result of the spline evaluation
        would be <= 1. We could set the ext class member of the spline to 3
        before reurning it, but event that would not be right in general.
        This is probably not a huge issue, but it should be addressed
        properly.
        """
        _x = self.xsection.x
        _y = numpy.exp(-column_density*self(_x))
        _fmt = dict(xname='Energy', xunits='keV', yname='Transmission factor')
        return xInterpolatedUnivariateSplineLinear(_x, _y, **_fmt)



def main():
    """Simple test code.
    """
    from ximpol.utils.matplotlib_ import pyplot as plt
    model = xpeInterstellarAbsorptionModel()
    trans = model.transmission_factor(1.e22)
    energy = numpy.linspace(1, 10, 10)
    print(trans(energy))
    plt.figure()
    model.xsection.plot(logx=True, logy=True, show=False)
    plt.figure()
    trans.plot(logx=True)
    

if __name__ == '__main__':
    main()

