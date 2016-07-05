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


from astropy.io import fits

from ximpol.utils.logging_ import logger
from ximpol.irf.base import OGIP_HEADER_SPECS
from ximpol.core.fitsio import xBinTableHDUBase
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear,\
    xInterpolatedBivariateSplineLinear


class xBinTableHDUSPECRESP(xBinTableHDUBase):

    """Binary table for the SPECRESP extension of a arf file.
    """

    NAME = 'SPECRESP'
    HEADER_KEYWORDS = [
        ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
        ('HDUCLAS2', 'SPECRESP', 'dataset contains spectral response')
    ] + OGIP_HEADER_SPECS
    DATA_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('SPECRESP', 'E', 'cm**2')
    ]


class xBinTableHDUVIGNETTING(xBinTableHDUBase):

    """Binary table for the VIGNETTING extension of a arf file.
    """

    NAME = 'VIGNETTING'
    HEADER_KEYWORDS = []
    DATA_SPECS = []

    def __init__(self, data, keywords=[], comments=[]):
        """Overloaded constructor.
        """
        energy, theta, vignetting = data
        ne = len(energy)
        nt = len(theta)
        assert vignetting.shape == (ne, nt)
        data = [energy.reshape((1, ne)),
                theta.reshape((1, nt)),
                vignetting.reshape((1, ne*nt))
                ]
        self.DATA_SPECS = [
            ('ENERGY'    , '%dE' % ne, 'keV'),
            ('THETA'     , '%dE' % nt, 'arcmin'),
            ('VIGNETTING', '%dE' % (ne*nt))
        ]
        self.HEADER_KEYWORDS = [
            ('TDIM3'     , '(%d, %d)' % (nt, ne))
        ]
        xBinTableHDUBase.__init__(self, data, keywords, comments)


class xEffectiveArea(xInterpolatedUnivariateSplineLinear):

    """Class describing an effectiva area.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    Arguments
    ---------
    arf_file_path : str
        The path to the .arf FITS file containing the effective area table.

    Example
    -------
    >>> import os
    >>> import numpy
    >>> from ximpol import XIMPOL_IRF
    >>>
    >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
    >>> aeff = xEffectiveArea(file_path)
    >>> x = numpy.arange(1, 10, 1)
    >>> print(aeff(x))
    >>> aeff.view()
    """

    def __init__(self, arf_file_path):
        """Constructor.
        """
        logger.info('Reading effective area data from %s...' % arf_file_path)
        self.hdu_list = fits.open(arf_file_path)
        self.hdu_list.info()
        _data = self.hdu_list['SPECRESP'].data
        _x = 0.5*(_data.field('ENERG_LO') + _data.field('ENERG_HI'))
        _y = _data.field('SPECRESP')
        fmt = dict(xname='Energy', xunits='keV', yname='Effective area',
                   yunits='cm$^2$')
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)
        _x = self.hdu_list['VIGNETTING'].data['ENERGY'][0]
        _y = self.hdu_list['VIGNETTING'].data['THETA'][0]
        _z = self.hdu_list['VIGNETTING'].data['VIGNETTING'][0]
        fmt = dict(xname='Energy', xunits='keV', yname='Off-axis angle',
                   yunits='arcmin', zname='Vignetting')
        self.vignetting = xInterpolatedBivariateSplineLinear(_x, _y, _z, **fmt)

    def eval_(self, energy, theta):
        """Return the effective area at a given energy and off-axis angle.
        """
        return self(energy)*self.vignetting(energy, theta)

    def view(self, off_axis_angle = 10., show=True):
        """Plot the effective area.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        plt.figure('Effective area')
        xInterpolatedUnivariateSplineLinear.plot(self, show=False,
                                                 label='On axis')
        plt.plot(self.x, self.eval_(self.x, off_axis_angle),
                 label='%s arcmin off-axis' % off_axis_angle)
        plt.legend(bbox_to_anchor=(0.85, 0.75))
        plt.figure('Vignetting')
        self.vignetting.plot(show=False)
        if show:
            plt.show()

            

if __name__ == '__main__':
    main()
