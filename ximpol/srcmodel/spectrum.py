#!/urs/bin/env python
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
from ximpol.core.rand import xUnivariateGenerator, xUnivariateAuxGenerator


class xSpectralComponent:

    def build_count_spectrum(self, effective_area, t):
        """
        """
        pass

    def build_light_curve(self):
        """
        """
        pass


class xCountSpectrum(xUnivariateAuxGenerator):

    """Class representing a count spectrum, i.e., the convolution of the
    source photon spectrum and the detector effective area.
    """

    def __init__(self, dNdE, aeff, t):
        """Constructor.
        """
        fmt = dict(auxname='Time', auxunits='s', rvname='Energy',
                   rvunits='keV',  pdfname='dN/dE',
                   pdfunits='cm$^{-2}$ s$^{-1}$ keV$^{-1}$')

        def _pdf(E, t):
            """Return the convolution between the effective area and
            the input photon spectrum.
            """
            return dNdE(E, t)*numpy.tile(aeff(E[0]), (t.shape[0], 1))

        xUnivariateAuxGenerator.__init__(self, t, aeff.x, _pdf, **fmt)
        self.light_curve = self.build_light_curve()

    def build_light_curve(self):
        """Build the light curve, i.e., a linear spline of the integral
        flux values as a function of time.
        """
        fmt = dict(rvname=self.xname, rvunits=self.xunits,
                   pdfname='Light curve', pdfunits='cm$^{-2}$ s$^{-1}$')
        _f = numpy.array([self.slice(_t).norm() for _t in self.x])
        return xUnivariateGenerator(self.x, _f, **fmt)





def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF
    from ximpol.irf.arf import xEffectiveArea

    def dNdE(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
    aeff = xEffectiveArea(file_path)
    t = numpy.linspace(0, 25, 100)
    c = xCountSpectrum(dNdE, aeff, t)
    c.light_curve.plot()
    c.plot()


if __name__ == '__main__':
    main()
