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


import numpy
from astropy.io import fits

from ximpol.utils.logging_ import logger
from ximpol.irf.base import OGIP_HEADER_SPECS
from ximpol.core.fitsio import xBinTableHDUBase
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.rand import xUnivariateAuxGenerator


class xBinTableHDUMATRIX(xBinTableHDUBase):

    """Binary table for the MATRIX extension of a rmf file.
    """

    NAME = 'MATRIX'
    HEADER_KEYWORDS = [
        ('HDUCLAS1', 'RESPONSE'  , 'dataset relates to spectral response'),
        ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
        ('CHANTYPE', 'PI '       , 'Detector Channel Type in use (PHA or PI)')
    ] + OGIP_HEADER_SPECS
    DATA_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('N_GRP'   , 'I'),
        ('F_CHAN'  , 'I'),
        ('N_CHAN'  , 'I'),
        ('MATRIX'  , None) #Mind this field is set at creation time
    ]

    def __init__(self, num_chans, data=None, keywords=[], comments=[]):
        """Overloaded constructor.
        """
        self.DATA_SPECS[-1] = ('MATRIX', '%dE' % num_chans)
        xBinTableHDUBase.__init__(self, data, keywords, comments)


class xBinTableHDUEBOUNDS(xBinTableHDUBase):

    """Binary table for the MATRIX extension of a rmf file.
    """

    NAME = 'EBOUNDS'
    HEADER_KEYWORDS = [
        ('CHANTYPE', 'PI'        , 'Channel type'),
        ('CONTENT' , 'Response Matrix', 'File content'),
        ('HDUCLAS1', 'RESPONSE'  , 'Extension contains response data  '),
        ('HDUCLAS2', 'EBOUNDS '  , 'Extension contains EBOUNDS')
    ] + OGIP_HEADER_SPECS
    DATA_SPECS = [
        ('CHANNEL', 'I'),
        ('E_MIN'  , 'E', 'keV'),
        ('E_MAX'  , 'E', 'keV')
    ]


class xEnergyDispersionMatrix(xUnivariateAuxGenerator):

    """Class encapsulating the energy dispersion matrix, as stored in the
    MATRIX extension of a .rmf file.

    In order to streamline performance, the energy grid is down-sampled
    to the value of the `num_aux_points` parameter, in such a way that
    the xUnivariateAuxGenerator.build_vppf()` doesn't take forever.

    Arguments
    ---------
    hdu : FITS hdu
       The MATRIX hdu in the .rmf FITS file.

    num_aux_point : int
       The number of points that the energy dispersion matrix should be\
       down-sampled to.

    Warning
    -------
    If the value of `num_aux_points` is too small, then the two-dimensional
    underlying spline tends to have blob-like features, and the vertical slices
    are no longer necessarily accurate representations of the energy dispersion.
    We should keep an eye on it.
    """

    def __init__(self, hdu, num_aux_points=200):
        """Constructor.
        """
        # First build a bivariate spline with the full data grid.
        _matrix = hdu.data
        _x = 0.5*(_matrix['ENERG_LO'] + _matrix['ENERG_HI'])
        _y = numpy.arange(0, len(_matrix['MATRIX'][0]), 1)
        _z = _matrix['MATRIX']
        _pdf = xInterpolatedBivariateSplineLinear(_y, _x, _z.transpose())
        # Then initialize the actual xUnivariateAuxGenerator object
        # with a down-sampled aux axis.
        _aux = numpy.linspace(_pdf.ymin(), _pdf.ymax(), num_aux_points)
        _rv = _y
        fmt = dict(auxname='Energy', auxunits='keV', rvname='Channel',
                   pdfname='Probability density')
        xUnivariateAuxGenerator.__init__(self, _aux, _rv, _pdf, **fmt)

    def rvs(self, aux):
        """Overloaded method.

        We want to return an integer, here.
        """
        val = xUnivariateAuxGenerator.rvs(self, aux)
        return numpy.ndarray.astype(numpy.rint(val), numpy.int16)


class xEnergyDispersionBounds(xInterpolatedUnivariateSplineLinear):

    """Class encapsulating the bounds for the energy dispersion matrix, as
    stored in the EBOUNDS extension of a .rmf file.
    """

    def __init__(self, hdu):
        """Constructor.
        """
        _bounds = hdu.data
        _x = _bounds['CHANNEL']
        _y = 0.5*(_bounds['E_MIN'] + _bounds['E_MAX'])
        fmt = dict(xname='Channel', yname='Energy', yunits='keV')
        xInterpolatedUnivariateSplineLinear.__init__(self, _x, _y, **fmt)


class xEnergyDispersion:

    """Class representing the energy dispersion.

    Arguments
    ---------
    rmf_file_path : str
        The path to the .rmf FITS file containing the energy dispersion tables.
    """

    def __init__(self, mrf_file_path):
        """Constructor.
        """
        logger.info('Reading energy dispersion data from %s...' % mrf_file_path)
        self.hdu_list = fits.open(mrf_file_path)
        self.hdu_list.info()
        self.matrix = xEnergyDispersionMatrix(self.hdu_list['MATRIX'])
        self.ebounds = xEnergyDispersionBounds(self.hdu_list['EBOUNDS'])

    def view(self, show=True):
        """Plot the energy dispersion.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        plt.figure('Energy redistribution')
        self.matrix.plot(show=False)
        plt.figure('Energy bounds')
        self.ebounds.plot(overlay=False, show=False)
        energy = 0.5*(self.matrix.xmin() + self.matrix.xmax())
        plt.figure('Energy dispersion @ %.2f keV' % energy)
        vslice = self.matrix.vslice(energy)
        vslice.plot(overlay=False, show=False)
        plt.text(0.1, 0.9, '$E = %.2f\\ \\rm{keV}$' % energy,
                 transform=plt.gca().transAxes)
        ppf = vslice.build_ppf()
        eres = 0.5*(ppf(0.8413) - ppf(0.1586))/ppf(0.5)
        plt.text(0.1, 0.85, '$\sigma_E/E = %.3f$' % eres,
                 transform=plt.gca().transAxes)
        if show:
            plt.show()



if __name__ == '__main__':
    main()
