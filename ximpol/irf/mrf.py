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
from ximpol.irf.base import xColDefsBase, OGIP_HEADER_SPECS
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


"""Header specifications for the MODFRESP extension of .mrf FITS files.
"""
MODFRESP_HEADER_SPECS = [
    ('EXTNAME' , 'MODFRESP', 'name of this binary table extension'),
    ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
    ('HDUCLAS2', 'MODFRESP', 'dataset contains modulation response')
] + OGIP_HEADER_SPECS


class xColDefsMODFRESP(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the MODFRESP extension
    of .mrf FITS files.
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('MODFRESP', 'E', None)
    ]


class xModulationFactor(xInterpolatedUnivariateSplineLinear):

    """Class describing the modulation factor.

    The effective area is essentially a linear spline, with built-in facilities
    for evaluation and plotting.

    Arguments
    ---------
    mrf_file_path : str
        The path to the .mrf FITS file containing the effective area table.

    Example
    -------
    >>> import os
    >>> import numpy
    >>> from ximpol import XIMPOL_IRF
    >>>
    >>> file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    >>> modf = xModulationFactor(file_path)
    >>> x = numpy.arange(1, 10, 1)
    >>> print(modf(x))
    >>> modf.plot(overlay=False)
    """

    def __init__(self, mrf_file_path):
        """Constructor.
        """
        logger.info('Reading modulation factor data from %s...' % mrf_file_path)
        hdu_list = fits.open(mrf_file_path)
        _data = hdu_list['MODFRESP'].data
        _x = 0.5*(_data.field('ENERG_LO') + _data.field('ENERG_HI'))
        _y = _data.field('MODFRESP')
        hdu_list.close()
        fmt = dict(xname='Energy', xunits='keV', yname='Modulation factor')
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)


def main():
    """
    """
    import os
    import numpy
    from ximpol import XIMPOL_IRF

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.mrf')
    modf = xModulationFactor(file_path)
    x = numpy.arange(1, 10, 1)
    print(modf(x))
    modf.plot(overlay=False)



if __name__ == '__main__':
    main()
