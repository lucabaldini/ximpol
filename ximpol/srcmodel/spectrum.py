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
from ximpol.utils.units_ import erg2keV
from ximpol.utils.logging_ import logger
from ximpol.srcmodel.gabs import xpeInterstellarAbsorptionModel


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


def pl_norm(integral, emin, emax, index, energy_power=0.):
    """Return the power-law normalization resulting in a given integral
    flux (or integral energy flux, or more in general integral of the
    flux multiplied by a generic power of the energy) between the minimum and
    maximum energies assuming a given spectral index.

    More specifically, given a power law of the form

    .. math::
       \\mathcal{S}(E) = C\\left( \\frac{E}{E_0} \\right)^{-\\Gamma}
       \\quad [{\\rm keV}^{-1}~{\\rm cm}^{-2}~{\\rm s}^{-1}],
    
    (where :math:`E_0 = 1~{\\rm keV}`) we define
    :math:`\\beta = (1 + p - \\Gamma)` and calculate
    
    .. math::
       I_{p} = \\int_{E_{\\rm min}}^{E_{\\rm max}} E^{p}\\mathcal{S}(E) dE =
       \\begin{cases}
       \\frac{C E_0^{\\Gamma}}{\\beta}
       \\left( E_{\\rm max}^{\\beta} - E_{\\rm min}^{\\beta}\\right)
       \\quad \\beta \\neq 0\\\\
       C E_0^{\\Gamma} \\ln \\left( E_{\\rm max}/E_{\\rm min} \\right)
       \\quad \\beta = 0\\\\
       \\end{cases}
       \\quad [{\\rm keV}^{p}~{\\rm cm}^{-2}~{\\rm s}^{-1}].
    
    Hence
    
    .. math::
        C =
        \\begin{cases}
        \\frac{I_p\\beta}{E_0^{\\Gamma}
        \\left( E_{\\rm max}^{\\beta} - E_{\\rm min}^{\\beta}\\right)}
        \\quad \\beta \\neq 0\\\\
        \\frac{I_p}{E_0^{\\Gamma}
        \\ln \\left( E_{\\rm max}/E_{\\rm min} \\right)}
        \\quad \\beta = 0.
        \\end{cases}

    Arguments
    ---------
    integral : float or array
        The value of the integral flux or energy-to-some-power flux

    emin : float
        The minimum energy for the integral flux

    emax : float
        The maximum energy for the integral flux

    index : float
        The power-law index

    energy_power : float
        The power of the energy in the integral

    erg : bool
        if True, convert erg to keV in the calculation.
    """
    assert emax > emin
    beta = 1 + energy_power - index
    if beta != 0:
        return integral*beta/(emax**beta - emin**beta)
    else:
        return integral/numpy.log(emax/emin)
    

def int_eflux2pl_norm(integral, emin, emax, index, erg=True):
    """Convert an integral energy flux into the corresponding power-law
    normalization.
    """
    if erg:
        integral = erg2keV(integral)
    return pl_norm(integral, emin, emax, index, 1.)


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

    def __init__(self, source_spectrum, aeff, t, column_density=0.,
                 redshift=0., scale_factor=1.):
        """Constructor.

        Arguments
        ---------
        source_spectrum : a Python function
            The source spectrum, i.e., a Python function that can be called
            with two numpy arrays E and t and returns the differential energy
            spectrum of the source at the appropriate enenrgy and time value(s).

        aeff : ximpol.irf.xEffectiveArea instance
            The instrument effective area.
        
        t : array
            The array of time values where we're sampling the light curve of the
            source.

        column_density : float, defaults to 0.
            The H column density to calculate the insterstellar absorption.

        redshift : float, defaults to 0.
            The source redshift.

        scale_factor : float, defaults to 1.
            An optional scale factor for the output count spectrum (see the
            warning below).

        Warning
        -------
        This was a workaround for running the MDP script on a pulsar,
        where we sample in phase and then have to multiply by the
        observation time.
        """

        def _pdf(E, t):
            """Return the convolution between the effective area and
            the input photon spectrum.

            Note that, due the way this callable is used within the constructor
            of xUnivariateAuxGenerator at this point E and t are two numpy
            arrays whose shape is `(num_time_samples, num_energy_samples)`.
            Particularly, `E[0]` is the array of energy values that we use for
            sampling the spectrum (nominally the same sampling that we use
            for the effective area).

            We start from all the stuff that is energy-independent, e.g., the
            effective area and the interstellar absorption and do the proper
            convolutions in energy space. Then we tile the values horizontally
            for as many times as the length of the array sampling the time,
            and at the end we multiply the resulting bidimensional array by
            an equal-size array containing the value of the differential
            energy spectrum on the time-energy grid.
            """
            # Grab the one-dimensional array used to sample the energy.
            _energy = E[0]
            # Calculate the effective area on the array.
            _vals = aeff(_energy)
            # Multiply by the scale factor.
            _vals *= scale_factor
            # If needed, fold in the transmission factor for the interstellar
            # absorption.
            if column_density > 0.:
                logger.info('Applying interstellar absorption (nH = %.3e)' %\
                            column_density)
                ism_model = xpeInterstellarAbsorptionModel()
                ism_trans = ism_model.transmission_factor(column_density)
                _vals *= ism_trans(_energy)
            # Tile the one dimensional vector horizontally.
            _vals = numpy.tile(_vals, (t.shape[0], 1))
            # And multiply by the source spectrum---we're done.
            _vals *= source_spectrum(E*(1 + redshift), t)
            return _vals

        fmt = dict(auxname='Time', auxunits='s', rvname='Energy',
                   rvunits='keV',  pdfname='dN/dE $\\times$ aeff',
                   pdfunits='s$^{-1}$ keV$^{-1}$')
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
        ebins = zip(energy_binning[:-1], energy_binning[1:])
        if len(energy_binning) > 2:
            ebins.append((energy_binning[0], energy_binning[-1]))
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
