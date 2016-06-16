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
from collections import OrderedDict

from ximpol.srcmodel.img import xFITSImage
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.evt.event import xMonteCarloEventList
from ximpol.core.spline import xInterpolatedUnivariateSplineLinear
from ximpol.utils.units_ import keV2erg, ergcms2mcrab
from ximpol.utils.logging_ import logger


DEFAULT_MAX_VALIDITY_TIME = 10000000.


class xModelComponentBase:

    """Base class for the source object.

    Note that the source identifier defaults to none and is typically assign
    after the fact when the source itself is added to a source model.

    Arguments
    ---------
    name : string
        The name of the source.

    energy_spectrum : function
        The function object representing the energy spectrum.

    polarization_degree : function
        The function object representing the polarization degree.

    polarization_angle : function
        The function object representing the polarization angle.

    column_density : float
        The value of the column density (in cm^-2) used to calculate the
        Galactic absorption. Defaults to 0.
    
    redshift : float
        The source redshift. Defaults to 0.

    min_validity_time : float
        The minimum validity time for the source model. Defaults to 0.

    max_validity_time : float
        The maximum validity time for the source model. Defaults to 10 Ms.

    identifier : int
        A unique identifier of the source within a ROI model. Defaults to None
        and is, generally, automatically assigned when building the ROI (i.e.,
        you don't have to worry about it, but it's handy to have in the
        output event list).
    """

    def __init__(self, name, energy_spectrum, polarization_degree,
                 polarization_angle, column_density=0., redshift=0.,
                 min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME, identifier=None):
        """Constructor.
        """
        self.name = name
        self.setup(energy_spectrum, polarization_degree, polarization_angle)
        self.column_density = column_density
        self.redshift = redshift
        self.min_validity_time = min_validity_time
        self.max_validity_time = max_validity_time
        self.identifier = identifier
        # We cache here the integral flux, since we don't want to (potentially)
        # recalculate the integral each time the __str__() method is called.
        emin = 2.
        emax = 8.
        t = self.min_validity_time
        self.flux_ergcms = self.integral_flux(emin, emax, t)
        self.flux_mcrab = ergcms2mcrab(self.flux_ergcms)
        self.flux_label = 'Flux @ t = %d: %.3e erg/cm2/s (%.2f mcrab)' %\
                          (t, self.flux_ergcms, self.flux_mcrab)

    def integral_flux(self, emin=2.0, emax=8.0, t=None, erg=True):
        """Return the integral source flux at a generic time.

        This is achieved by taking a "slice" of the source spectrum
        at that time and integrating between a minimum and maximum energy.

        Arguments
        ---------
        emin : float
            The minimum integration energy (default 2 keV).

        emax : float
            The maximum integration energy (default 2 keV).

        t : float
            The time (default is the minimum source validity time).
        """
        if t is None:
            t = self.min_validity_time
        _x = numpy.linspace(emin, emax, 100)
        _y = _x*self.energy_spectrum(_x, t)
        _flux = xInterpolatedUnivariateSplineLinear(_x, _y).integral(emin, emax)
        if erg:
            _flux = keV2erg(_flux)
        return _flux

    @classmethod
    def sampling_time(cls, tstart, tstop):
        """Return the sampling time used for the lightcurve when generating
        events.

        Warning
        -------
        We have an issue open for this and at some point we should go back to
        it and make sure it is properly addressed.
        """
        return numpy.linspace(tstart, tstop, 100)

    def set_energy_spectrum(self, energy_spectrum):
        """Set the energy spectrum for the model component.

        Arguments
        ---------
        energy_spectrum : function
            The function object representing the energy spectrum.
        """
        self.energy_spectrum = energy_spectrum

    def set_polarization_degree(self, polarization_degree):
        """Set the polarization degree for the model component.

        Arguments
        ---------
        polarization_degree : function
            The function object representing the polarization degree.
        """
        self.polarization_degree = polarization_degree

    def set_polarization_angle(self, polarization_angle):
        """Set the polarization angle for the model component.

        Arguments
        ---------
        polarization_angle : function
            The function object representing the polarization angle.
        """
        self.polarization_angle = polarization_angle

    def setup(self, energy_spectrum, polarization_degree, polarization_angle):
        """Setup the model component in terms of energy spectrum and
        polarization degree and angle.
        """
        self.set_energy_spectrum(energy_spectrum)
        self.set_polarization_degree(polarization_degree)
        self.set_polarization_angle(polarization_angle)

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is a do-nothing function and should be re-implemented by
        each derived class.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        pass

    def __str__(self):
        """String formatting.

        Note that this base class has no ra and dec attributes, since they're
        not really relevant for extended sources initiazed from FITS images,
        but since all the other derived classed do have such attributes,
        check whether that's the case and, if necessary, print out the
        source position. This avoids code duplication in the derived classes.
        """
        text = '%s "%s" (id = %s)'%\
               (self.__class__.__name__, self.name, self.identifier)
        text += '\n    Validity time: [%f--%f]'  %\
                (self.min_validity_time, self.max_validity_time)
        text += '\n    Galactic column density: %.3e cm^{-2}' %\
                self.column_density
        text += '\n    Redshift: %.3f' % self.redshift
        text += '\n    %s' % self.flux_label
        if hasattr(self, 'ra') and hasattr(self, 'dec'):
            text += '\n    Position: RA = %s deg, Dec = %s deg' %\
                    (self.ra, self.dec)
        return text

    def rvs_event_list(self, aeff, psf, modf, edisp, **kwargs):
        """Extract a random event list for the model component.
        """
        # Create the event list and the count spectrum.
        event_list = xMonteCarloEventList()
        tsamples = self.sampling_time(kwargs['tstart'], kwargs['tstop'])
        logger.info('Sampling times: %s' % tsamples)
        count_spectrum = xCountSpectrum(self.energy_spectrum, aeff, tsamples,
                                        self.column_density, self.redshift)
        # Extract the number of events to be generated based on the integral
        # of the light curve over the simulation time.
        num_events = numpy.random.poisson(count_spectrum.light_curve.norm())
        logger.info('About to generate %d events...' % num_events)
        # Extract the event times and sort them.
        col_time = count_spectrum.light_curve.rvs(num_events)
        col_time.sort()
        event_list.set_column('TIME', col_time)
        # Extract the MC energies and smear them with the energy dispersion.
        col_mc_energy = count_spectrum.rvs(col_time)
        event_list.set_column('MC_ENERGY', col_mc_energy)
        col_pha = edisp.matrix.rvs(col_mc_energy)
        event_list.set_column('PHA', col_pha)
        event_list.set_column('ENERGY', edisp.ebounds(col_pha))
        # Extract the MC sky positions and smear them with the PSF.
        col_mc_ra, col_mc_dec = self.rvs_sky_coordinates(num_events)
        event_list.set_column('MC_RA', col_mc_ra)
        event_list.set_column('MC_DEC', col_mc_dec)
        col_ra, col_dec = psf.smear(col_mc_ra, col_mc_dec)
        event_list.set_column('RA', col_ra)
        event_list.set_column('DEC', col_dec)
        # Extract the photoelectron emission directions.
        pol_degree = self.polarization_degree(col_mc_energy, col_time,
                                              col_mc_ra, col_mc_dec)
        pol_angle = self.polarization_angle(col_mc_energy, col_time,
                                              col_mc_ra, col_mc_dec)
        col_pe_angle = modf.rvs_phi(col_mc_energy, pol_degree, pol_angle)
        event_list.set_column('PE_ANGLE', col_pe_angle)
        # Set the source ID.
        event_list.set_column('MC_SRC_ID', self.identifier)
        # Set the phase to rnd [0-1] for all non-periodic sources.
        phase=numpy.random.uniform(0,1,len(col_pe_angle))
        event_list.set_column('PHASE', phase)
        return event_list


class xPointSource(xModelComponentBase):

    """Class representing a steady point source.
    See :py:class:`ximpol.srcmodel.roi.xModelComponentBase` for the signature
    of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the source (in decimal degrees).

    dec : float
        The declination of the source (in decimal degrees).
    """

    def __init__(self, name, ra, dec, energy_spectrum, polarization_degree,
                 polarization_angle, column_density=0., redshift=0.,
                 min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME, identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, energy_spectrum,
                                     polarization_degree, polarization_angle,
                                     column_density, redshift,
                                     min_validity_time, max_validity_time,
                                     identifier)
        self.ra = ra
        self.dec = dec

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        ra = numpy.full(size, self.ra)
        dec = numpy.full(size, self.dec)
        return (ra, dec)


class xEphemeris:

    """Convenience class encapsulating a pulsar ephemeris.
    """

    def __init__(self, t0, nu0, nudot=0., nuddot=0., min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME):
        """Constructor.
        """
        self.t0 = t0
        self.nu0 = nu0
        self.nudot = nudot
        self.nuddot = nuddot
        self.min_validity_time = min_validity_time
        self.max_validity_time = max_validity_time

    def nu(self, t):
        """Return the source frequency at a given time.
        """
        assert (t >= self.min_validity_time) and (t <= self.max_validity_time)
        dt = (t - self.t0)
        return self.nu0 + self.nudot*dt + 0.5*self.nuddot*numpy.power(dt, 2.)

    def period(self, t):
        """Return the source period at a given time.
        """
        return 1./self.nu(t)

    def __str__(self):
        """String formatting.
        """
        return 't0 = %s s, nu0 = %s Hz, nudot = %s Hz/s, nuddot = %s Hz/s^2' %\
            (self.t0, self.nu0, self.nudot, self.nuddot)


class xPeriodicPointSource(xPointSource):

    """Class representing a periodic point source (e.g., a pulsar).
    See :py:class:`ximpol.srcmodel.roi.xModelComponentBase` for the signature
    of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the source (in decimal degrees).

    dec : float
        The declination of the source (in decimal degrees).

    ephemeris : :py:class:`ximpol.srcmodel.roi.xEphemeris object
        The source ephemeris.
    """

    def __init__(self, name, ra, dec, energy_spectrum, polarization_degree,
                 polarization_angle, ephemeris, column_density=0., redshift=0.,
                 min_validity_time=None, max_validity_time=None,
                 identifier=None):
        """Constructor.
        """
        if min_validity_time is None:
            min_validity_time = ephemeris.min_validity_time
        if max_validity_time is None:
            max_validity_time = ephemeris.max_validity_time
        xPointSource.__init__(self, name, ra, dec, energy_spectrum,
                              polarization_degree, polarization_angle,
                              column_density, redshift, min_validity_time,
                              max_validity_time, identifier)
        self.ephemeris = ephemeris

    def rvs_event_list(self, aeff, psf, modf, edisp, **kwargs):
        """Extract a random event list for the model component.

        TODO: here we should pass the sampling phase, instead?

        TODO: properly take into account the derivatives in the ephemeris.
        """
        # Create the event list and the count spectrum.
        event_list = xMonteCarloEventList()
        # Mind the count spectrum is made in phase!
        sampling_phase = numpy.linspace(0., 1., 100)
        count_spectrum = xCountSpectrum(self.energy_spectrum, aeff,
                                        sampling_phase, self.column_density,
                                        self.redshift)
        # All this is not properly taking into account the ephemeris.
        min_time = kwargs['tstart']
        max_time = kwargs['tstop']
        #min_time=sampling_time[0]
        #max_time = sampling_time[-1]
        delta_time = (max_time - min_time)
        period = self.ephemeris.period(min_time)
        # This is not accurate, as we are effectively discarding the last
        # fractional period. Need to think about it.
        num_periods = int(delta_time/period)
        num_expected_events = delta_time*count_spectrum.light_curve.norm()
        # Extract the number of events to be generated based on the integral
        # of the light curve over the simulation time.
        num_events = numpy.random.poisson(num_expected_events)
        # Extract the event phases and sort them.
        col_phase = count_spectrum.light_curve.rvs(num_events)
        event_list.set_column('PHASE', col_phase)
        col_period = numpy.random.randint(0, num_periods, num_events)
        col_time = (col_period + col_phase)*period
        event_list.set_column('TIME', col_time)
        # Extract the MC energies and smear them with the energy dispersion.
        col_mc_energy = count_spectrum.rvs(col_phase)
        event_list.set_column('MC_ENERGY', col_mc_energy)
        col_pha = edisp.matrix.rvs(col_mc_energy)
        event_list.set_column('PHA', col_pha)
        event_list.set_column('ENERGY', edisp.ebounds(col_pha))
        # Extract the MC sky positions and smear them with the PSF.
        col_mc_ra, col_mc_dec = self.rvs_sky_coordinates(num_events)
        event_list.set_column('MC_RA', col_mc_ra)
        event_list.set_column('MC_DEC', col_mc_dec)
        col_ra, col_dec = psf.smear(col_mc_ra, col_mc_dec)
        event_list.set_column('RA', col_ra)
        event_list.set_column('DEC', col_dec)
        # Extract the photoelectron emission directions.
        pol_degree = self.polarization_degree(col_mc_energy, col_phase,
                                              col_mc_ra, col_mc_dec)
        pol_angle = self.polarization_angle(col_mc_energy, col_phase,
                                            col_mc_ra, col_mc_dec)
        col_pe_angle = modf.rvs_phi(col_mc_energy, pol_degree, pol_angle)
        event_list.set_column('PE_ANGLE', col_pe_angle)
        # Set the source ID.
        event_list.set_column('MC_SRC_ID', self.identifier)
        event_list.sort()
        return event_list

    def __str__(self):
        """String formatting.
        """
        text = xPointSource.__str__(self)
        text += '\n    Ephemeris: %s' % self.ephemeris
        return text


class xUniformDisk(xModelComponentBase):

    """Class representing a uniform disk.
    See :py:class:`ximpol.srcmodel.roi.xModelComponentBase` for the signature
    of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the disk center (in decimal degrees).

    dec : float
        The declination of the disk center (in decimal degrees).

    radius : float
        The radius of the disk (in degrees).
    """

    def __init__(self, name, ra, dec, radius, energy_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME, identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, energy_spectrum,
                                     polarization_degree, polarization_angle,
                                     column_density, redshift,
                                     min_validity_time, max_validity_time,
                                     identifier)
        self.ra = ra
        self.dec = dec
        self.radius = radius

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        The algorithm is taken from
        http://mathworld.wolfram.com/DiskPointPicking.html

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        r = self.radius*numpy.sqrt(numpy.random.sample(size))
        theta = numpy.random.uniform(0, 2*numpy.pi, size)
        ra = self.ra + r*numpy.cos(theta)
        dec = self.dec + r*numpy.sin(theta)
        return (ra, dec)

    def __str__(self):
        """String formatting.
        """
        text = xModelComponentBase.__str__(self)
        text += '\n    Radius: %s deg' % self.radius
        return text


class xGaussianDisk(xModelComponentBase):

    """Class representing a (azimuthally simmetric) gaussian disk.
    See :py:class:`ximpol.srcmodel.roi.xModelComponentBase` for the signature
    of the base class.

    Arguments
    ---------
    ra : float
        The right ascension of the disk center (in decimal degrees).

    dec : float
        The declination of the disk center (in decimal degrees).

    sigma : float
        The root mean square of the disk (in degrees).
    """

    def __init__(self, name, ra, dec, sigma, energy_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME, identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, energy_spectrum,
                                     polarization_degree, polarization_angle,
                                     column_density, redshift,
                                     min_validity_time, max_validity_time,
                                     identifier)
        self.ra = ra
        self.dec = dec
        self.sigma = sigma
        self.__mean = [self.ra, self.dec]
        self.__cov = [[sigma**2., 0.], [0., sigma**2.]]

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        This is returning an array of the proper length with identical values.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        rvs = numpy.random.multivariate_normal(self.__mean, self.__cov, size)
        ra, dec = rvs[:,0], rvs[:,1]
        return (ra, dec)

    def __str__(self):
        """String formatting.
        """
        text = xModelComponentBase.__str__(self)
        text += '\n    Sigma: %s deg' % self.sigma
        return text


class xExtendedSource(xModelComponentBase):

    """Class representing an extended source.
    See :py:class:`ximpol.srcmodel.roi.xModelComponentBase` for the signature
    of the base class.

    Arguments
    ---------
    img_file_path : string
        The path to the FITS file containing the image of the source.
    """

    def __init__(self, name, img_file_path, energy_spectrum,
                 polarization_degree, polarization_angle, column_density=0.,
                 redshift=0., min_validity_time=0.,
                 max_validity_time=DEFAULT_MAX_VALIDITY_TIME, identifier=None):
        """Constructor.
        """
        xModelComponentBase.__init__(self, name, energy_spectrum,
                                     polarization_degree, polarization_angle,
                                     column_density, redshift,
                                     min_validity_time, max_validity_time,
                                     identifier)
        self.image = xFITSImage(img_file_path)

    def rvs_sky_coordinates(self, size=1):
        """Generate random coordinates for the model component.

        Arguments
        ---------
        size : float
            The number of sky coordinate pairs to be generated.
        """
        return self.image.rvs_coordinates(size)

    def __str__(self):
        """String formatting.
        """
        text = xModelComponentBase.__str__(self)
        text += '\n    FITS image: %s' % self.image
        return text


class xROIModel(OrderedDict):

    """Class describing a full ROI (region of interest) model.

    This is essentially an (ordered) collection of component objects
    (i.e., instances of classes inheriting from xModelComponentBase)
    than can be accessed by source name.

    Arguments
    ---------
    ra_center : float
        The right ascension of the center of the ROI (in decimal degrees).

    dec_center : float
        The declination of the center of the ROI (in decimal degrees).
    """

    def __init__(self, ra_center, dec_center):
        """Constructor.
        """
        OrderedDict.__init__(self)
        self.ra = ra_center
        self.dec = dec_center

    def add_source(self, source):
        """Add a source to the ROI.
        """
        source.identifier = len(self)
        self[source.name] = source

    def add_sources(self, *sources):
        """Add an arbitrary number of sources to the ROI.
        """
        for source in sources:
            self.add_source(source)

    def __add__(self, other):
        """Combine different ROI models.
        """
        assert self.__class__.__name__ == other.__class__.__name__
        roi_model = xROIModel(self.ra, self.dec)
        for source in self.values():
            roi_model.add_source(source)
        for source in other.values():
            roi_model.add_source(source)
        return roi_model

    def min_validity_time(self):
        """Return the minimum validity time for the ROI model.
        """
        return max([source.min_validity_time for source in self.values()])

    def max_validity_time(self):
        """Return the maximum validity time for the ROI model.
        """
        return min([source.max_validity_time for source in self.values()])

    def __str__(self):
        """String formatting.
        """
        txt = 'ROI centered at (%.4f, %.4f):\n' % (self.ra, self.dec)
        for source in self.values():
            txt += '- %s\n' % source
        return txt.strip('\n')

    def build_hdu(self):
        """Build a FITS HDU for the source model.

        This can be used to be written in the output file.
        """
        pass

    def rvs_event_list(self, aeff, psf, modf, edisp, **kwargs):
        """Extract an event list for the full ROI.

        Arguments
        ---------
        aeff : :py:class:`ximpol.irf.arf.xEffectiveArea` object.
            The effective area to be used.

        psf : :py:class:`ximpol.irf.psf.xPointSpreadFunction` object.
            The PSF to be used.

        modf : :py:class:`ximpol.irf.mrf.xModulationFactor` object.
            The modulation factor to the used.

        edisp : :py:class:`ximpol.irf.rmf.xEnergyDispersion` object.
            The energy dispersion to be used.

        sampling_time : array
            The array to sample the source light curve.

        Warning
        -------
        The sampling_time should not be the same for all sources, and each
        source should be able to decide its own in a sensible way.
        (See issue #44.)
        """
        event_list = xMonteCarloEventList()
        for source in self.values():
            event_list += source.rvs_event_list(aeff, psf, modf, edisp,
                                                **kwargs)
        event_list.sort()
        return event_list



def main():
    """Simple test code.
    """
    import os
    from ximpol.srcmodel.spectrum import power_law
    from ximpol.srcmodel.polarization import constant
    from ximpol import XIMPOL_CONFIG
    src_base_model = dict(energy_spectrum = power_law(1., 2.),
                          polarization_degree = constant(0.1),
                          polarization_angle = constant(20.),
                          column_density = 1.0e-22,
                          redshift = 0.0012,
                          identifier = 0)
    s1 = xModelComponentBase('Source 1', **src_base_model)
    print(s1)
    s2 = xPointSource('Source 2', ra=10., dec=15., **src_base_model)
    print(s2)
    _ephem = xEphemeris(10., 0.033, min_validity_time=10.,
                        max_validity_time=10000.)
    s3 = xPeriodicPointSource('Source 3', ra=10., dec=15., ephemeris=_ephem,
                              **src_base_model)
    print(s3)
    s4 = xUniformDisk('Source 4', ra=10., dec=15., radius=0.01,
                      **src_base_model)
    print(s4)
    s5 = xGaussianDisk('Source 5', ra=10., dec=15., sigma=0.01,
                       **src_base_model)
    print(s5)
    img_file_path = os.path.join(XIMPOL_CONFIG, 'fits',
                                 'crab_0p3_10p0_keV.fits')
    s6 = xExtendedSource('Source 6', img_file_path, **src_base_model)
    print(s6)
    

if __name__ == '__main__':
    main()
