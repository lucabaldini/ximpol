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
from ximpol.irf.base import xColDefsBase, OGIP_HEADER_SPECS
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear


"""Header specifications for the MATRIX extension of .rmf FITS files.
"""
MATRIX_HEADER_SPECS = [
    ('EXTNAME' , 'MATRIX'    , 'name of this binary table extension'),
    ('HDUCLAS1', 'RESPONSE'  , 'dataset relates to spectral response'),
    ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
    ('DETCHANS', '1024'      , 'total number of detector channels'),
    ('CHANTYPE', 'PI '       , 'Detector Channel Type in use (PHA or PI)')
] + OGIP_HEADER_SPECS


"""Header specifications for the EBOUNDS extension of .rmf FITS files.
"""
EBOUNDS_HEADER_SPECS = [
    ('EXTNAME' , 'EBOUNDS '  , 'Name of this binary table extension'),
    ('CHANTYPE', 'PI'        , 'Channel type'),
    ('CONTENT' , 'Response Matrix', 'File content'),
    ('HDUCLAS1', 'RESPONSE'  , 'Extension contains response data  '),
    ('HDUCLAS2', 'EBOUNDS '  , 'Extension contains EBOUNDS'),
    ('DETCHANS', 1024        , 'Total number of detector channels')
] + OGIP_HEADER_SPECS



class xColDefsMATRIX(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the MATRIX extension
    of .rmf FITS files.
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('N_GRP'   , 'I', None),
        ('F_CHAN'  , 'I', None),
        ('N_CHAN'  , 'I', None),
        ('MATRIX'  , '1024E', None)
    ]


class xColDefsEBOUNDS(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the EBOUNDS extension
    of .rmf FITS files.
    """

    COLUMN_SPECS = [
        ('CHANNEL', 'I', None),
        ('E_MIN'  , 'E', 'keV'),
        ('E_MAX'  , 'E', 'keV')
    ]


class xEnergyDispersionMatrix(xInterpolatedBivariateSplineLinear):

    """Class encapsulating the energy dispersion matrix, as stored in the
    MATRIX extension of a .rmf file.

    This is essentially a bivariate linear spline on a rectangular mesh.
    """

    def __init__(self, hdu):
        """Constructor.
        """
        _matrix = hdu.data
        _x = 0.5*(_matrix['ENERG_LO'] + _matrix['ENERG_HI'])
        _y = numpy.arange(0, len(_matrix['MATRIX'][0]), 1)
        _z = _matrix['MATRIX']
        fmt = dict(xname='Energy', xunits='keV', yname='Channel',
                   zname='Probability density')
        xInterpolatedBivariateSplineLinear.__init__(self, _x, _y, _z, **fmt)


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
        The path to the .mrf FITS file containing the energy dispersion tables.
    """

    def __init__(self, mrf_file_path):
        """Constructor.
        """
        logger.info('Reading energy dispersion data from %s...' % mrf_file_path)
        hdu_list = fits.open(mrf_file_path)
        hdu_list.info()
        self.matrix = xEnergyDispersionMatrix(hdu_list['MATRIX'])
        self.ebounds = xEnergyDispersionBounds(hdu_list['EBOUNDS'])
        hdu_list.close()

    def plot(self):
        """Plot the energy dispersion.
        """
        from ximpol.utils.matplotlib_ import pyplot as plt
        from ximpol.utils.matplotlib_ import context_two_by_two
        emin = self.matrix.xmin()
        emax = self.matrix.xmax()

        def _plot_vslice(energy, position):
            """Convenience function to plot a generic vertical slice of the
            energy dispersion.
            """
            ax = plt.subplot(2, 2, position)
            vslice = self.matrix.vslice(energy)
            vslice.plot(overlay=False, show=False)
            plt.text(0.1, 0.9, '$E = %.2f\\ \\rm{keV}$' % energy,
                     transform=ax.transAxes)
            ppf = vslice.build_ppf()
            eres = 0.5*(ppf(0.8413) - ppf(0.1586))/ppf(0.5)
            plt.text(0.1, 0.85, '$\sigma_E/E = %.3f$' % eres,
                     transform=ax.transAxes)

        with context_two_by_two():
            plt.figure(1)
            ax = plt.subplot(2, 2, 1)
            self.matrix.plot(show=False)
            ax = plt.subplot(2, 2, 2)
            self.ebounds.plot(overlay=False, show=False)
            _plot_vslice(emin + 0.333*(emax - emin), 3)
            _plot_vslice(emin + 0.666*(emax - emin), 4)
        plt.show()


def main():
    """
    """
    import os
    import numpy
    from ximpol import XIMPOL_IRF

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.rmf')
    edisp = xEnergyDispersion(file_path)
    edisp.plot()


if __name__ == '__main__':
    main()
