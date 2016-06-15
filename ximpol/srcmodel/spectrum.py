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


import numpy
from ximpol.core.rand import xUnivariateGenerator, xUnivariateAuxGenerator
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.core.spline import xInterpolatedBivariateSplineLinear
from ximpol.irf.mrf import mdp99, xMDPTable


def power_law(C, Gamma):
    """Photon energy spectrum as a function of energy and time.

    If C and Gamma are callable, we assume that the argument of the __call__
    function is the time, and this is how we treat them internally.
    """
    if hasattr(C, '__call__') and hasattr(Gamma, '__call__'):
        def _function(E, t):
            return C(t)*numpy.power(E, -Gamma(t))
    else:
        def _function(E, t):
            return C*numpy.power(E, -Gamma)
    return _function


class xCountSpectrum(xUnivariateAuxGenerator):

    """Class representing a count spectrum, i.e., the convolution of the
    source photon spectrum and the detector effective area

    .. math::
       \\mathcal{C}(E, t) = \\mathcal{S}(E, t) \\times A_{\\rm eff}(E)
       \\quad [\\text{s}^{-1}~\\text{keV}^{-1}].

    This is a subclass of xUnivariateAuxGenerator, providing all the facilities
    implemented in a bivariate spline, along with the capability of
    extracting random numbers.

    Note that the light curve corresponding to the count spectrum is
    calculated when a class object is instantiated.
    """

    def __init__(self, source_spectrum, aeff, t, scale=1.):
        """Constructor.

        Warning
        -------
        Do we really want the option to pass a scale here?
        This was a workaround for running the MDP script on a pulsar,
        where we sample in phase and then have to multiply by the
        number of periods.
        """
        fmt = dict(auxname='Time', auxunits='s', rvname='Energy',
                   rvunits='keV',  pdfname='dN/dE $\\times$ aeff',
                   pdfunits='s$^{-1}$ keV$^{-1}$')

        def _pdf(E, t):
            """Return the convolution between the effective area and
            the input photon spectrum.
            """
            return scale*source_spectrum(E, t)*numpy.tile(aeff(E[0]),
                                                          (t.shape[0], 1))

        xUnivariateAuxGenerator.__init__(self, t, aeff.x, _pdf, **fmt)
        self.light_curve = self.build_light_curve()

    def build_time_integral(self, tmin=None, tmax=None):
        """Build the time-integrated count spectrum, i.e.

        .. math::
           \\int_{t_{\\rm min}}^{t_{\\rm max}} \\mathcal{C}(E, t) dt \\quad
           [\\text{keV}^{-1}].

        The output is stored in the form of a xUnivariateGenerator, featuring
        all the spline facilities, along with the capability of extracting
        random numbers.
        """
        if tmin is None:
            tmin = self.xmin()
        if tmax is None:
            tmax = self.xmax()
        _x = self.y
        _y = numpy.array([self.hslice(_E).integral(tmin, tmax) for _E in _x])
        fmt = dict(rvname=self.yname, rvunits=self.yunits,
                   pdfname='Time-integrated (%d--%d s) spectrum' %\
                   (tmin, tmax), pdfunits='keV$^{-1}$')
        return xUnivariateGenerator(_x, _y, **fmt)

    def build_energy_integral(self, emin=None, emax=None):
        """Build the energy-integrated count spectrum, i.e.

        .. math::
           \\int_{E_{\\rm min}}^{E_{\\rm max}} \\mathcal{C}(E, t) dE \\quad
           [\\text{Hz}].

        The output is stored in the form of a xUnivariateGenerator, featuring
        all the spline facilities, along with the capability of extracting
        random numbers.
        """
        if emin is None:
            emin = self.ymin()
        if emax is None:
            emax = self.ymax()
        _x = self.x
        _y = numpy.array([self.vslice(_t).integral(emin, emax) for _t in _x])
        fmt = dict(rvname=self.xname, rvunits=self.xunits,
                   pdfname='Energy-integrated (%.2f--%.2f keV) spectrum' %\
                   (emin, emax), pdfunits='Hz')
        return xUnivariateGenerator(_x, _y, **fmt)

    def build_light_curve(self):
        """Build the light curve, i.e., the count spectrum, integrated over the
        entire energy range, as a function of time.
        """
        return self.build_energy_integral(self.ymin(), self.ymax())

    def num_expected_counts(self, tmin=None, tmax=None, emin=None, emax=None):
        """Return the number of expected counts within a given time interval
        and energy range

        .. math::
           \\int_{t_{\\rm min}}^{t_{\\rm max}}
           \\int_{E_{\\rm min}}^{E_{\\rm max}}
           \\mathcal{C}(E, t) dt dE

        """
        if tmin is None:
            tmin = self.xmin()
        if tmax is None:
            tmax = self.xmax()
        if emin is None:
            emin = self.ymin()
        if emax is None:
            emax = self.ymax()
        return self.integral(tmin, tmax, emin, emax)

    def build_mdp_table(self, energy_binning, modulation_factor):
        """Calculate the MDP values in energy bins, given the modulation
        factor of the instrument as a function of the energy.

        Arguments
        ---------
        energy_binning : array
            The energy binning
        modulation_factor : ximpol.irf.mrf.xModulationFactor instance
            The instrument modulation factor as a function of the energy.
        """
        # Build the time-integrated spectrum.
        time_integrated_spectrum = self.build_time_integral()
        # Build the modulation-factor spectrum, i.e. the object that we
        # integrate to calculate the effective modulation factor in a
        # given energy bin for a given spectrum.
        _x = time_integrated_spectrum.x
        _y = time_integrated_spectrum.y*modulation_factor(_x)
        mu_spectrum = xInterpolatedUnivariateSplineLinear(_x, _y)
        # Loop over the energy bins and calculate the actual MDP values.
        # Note that we also calculate the MDP on the entire energy range.
        observation_time = self.xmax() - self.xmin()
        mdp_table = xMDPTable(observation_time)
        ebins = zip(energy_binning[:-1], energy_binning[1:]) +\
                [(energy_binning[0], energy_binning[-1])]
        for _emin, _emax in ebins:
            num_counts = self.num_expected_counts(emin=_emin, emax=_emax)
            mu_eff = mu_spectrum.integral(_emin, _emax)/num_counts
            mdp = mdp99(mu_eff, num_counts)
            mdp_table.add_row(mdp, _emin, _emax, mu_eff, num_counts)
        return mdp_table


def main():
    """
    """
    import os
    from ximpol import XIMPOL_IRF
    from ximpol.irf.arf import xEffectiveArea

    def source_spectrum(E, t):
        """Function defining a time-dependent energy spectrum.
        """
        return 10.0*(1.0 + numpy.cos(t))*numpy.power(E, (-2.0 + 0.01*t))

    file_path = os.path.join(XIMPOL_IRF,'fits','xipe_baseline.arf')
    aeff = xEffectiveArea(file_path)
    t = numpy.linspace(0, 25, 100)
    c = xCountSpectrum(source_spectrum, aeff, t)
    c.light_curve.plot()
    c.plot()


if __name__ == '__main__':
    main()
