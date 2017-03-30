#!/usr/bin/env python
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


"""Module encapsulating the FITS spectra structure and related facilities.
"""

import numpy
from astropy.io import fits
from astropy import wcs

from ximpol.utils.logging_ import logger, abort
from ximpol.evt.event import xEventFile
from ximpol.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ximpol.irf.mrf import xAzimuthalResponseGenerator
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.srcmodel.img import xFITSImage
from ximpol import xpColor
from ximpol.irf import irf_file_path



class xEventBinningBase:

    """Base class for the event binning.
    """

    def __init__(self, file_path, **kwargs):
        """Constructor.
        """
        self.event_file = xEventFile(file_path)
        self.event_data = self.event_file.event_data
        self.kwargs = kwargs
        self.process_kwargs()

    def get(self, key, default=None):
        """Convenience method to address the keyword aguments.
        """
        return self.kwargs.get(key, default)

    def set(self, key, value):
        """Convenience method to set keyword arguments.
        """
        logger.info('Setting %s to %s...' % (key, value))
        self.kwargs[key] = value

    def build_primary_hdu(self):
        """Build the primary HDU for the output file.
        """
        primary_hdu = xPrimaryHDU()
        primary_hdu.setup_header(self.event_file.primary_keywords())
        primary_hdu.add_keyword('BINALG', self.get('algorithm'),
                                'the binning algorithm used')
        primary_hdu.add_comment('%s run with kwargs %s' %\
                                (self.__class__.__name__, self.kwargs))
        return primary_hdu

    @classmethod
    def read_binning(cls, file_path):
        """Read a custom binning from file and return a numpy array.
        """
        return numpy.loadtxt(file_path)

    @classmethod
    def bin_centers(cls, bin_edges):
        """Return an array of bin centers given an array of bin edges.

        Arguments
        ---------
        bin_edges : 1-d array of length (n + 1).
            The array with the bin edges.

        Returns
        -------
        1-d array of length n.
            The array with the values of the bin centers.
        """
        assert bin_edges.ndim == 1
        return 0.5*(bin_edges[:-1] + bin_edges[1:])

    @classmethod
    def bin_widths(cls, bin_edges):
        """Return an array of bin widths given an array of bin edges.

        Arguments
        ---------
        bin_edges : 1-d array of length (n + 1).
            The array with the bin edges.

        Returns
        -------
        1-d array of length n.
            The array with the values of the bin widths.
        """
        assert bin_edges.ndim == 1
        return (bin_edges[1:] - bin_edges[:-1])

    @classmethod
    def bin_edge_pairs(cls, bin_edges, uber=False):
        """Return a list of (min, max) tuples of length 2 given an array of
        bin edges. This is handy when one needs to loop explicitly over the
        intervals defined by the array bounds, e.g., when running xpselect
        in specific energy, time or phase intervals.

        Arguments
        ---------
        bin_edges : 1-d array of length (n + 1).
            The array with the bin edges.

        uber : bool
            If true, add the overall interval (i.e., a tuple with the absolute
            minimum and maximum in the binning) at the end.

        Example
        -------
        >>> import numpy
        >>> from ximpol.evt.binning import xEventBinningBase
        >>>
        >>> energy_binning = numpy.array([2., 4., 8.])
        >>> print(xEventBinningBase.bin_edge_pairs(energy_binning))
        >>> [(2.0, 4.0), (4.0, 8.0)]
        >>> print(xEventBinningBase.bin_edge_pairs(energy_binning, uber=True))
        >>> [(2.0, 4.0), (4.0, 8.0), (2.0, 8.0)]

        """
        assert bin_edges.ndim == 1
        pairs = zip(bin_edges[:-1], bin_edges[1:])
        if uber and len(pairs) > 1:
            pairs.append((bin_edges[0], bin_edges[-1]))
        return pairs

    @classmethod
    def equipopulated_binning(cls, num_bins, vector, min_value=None,
                              max_value=None):
        """
        """
        if min_value is None:
            min_value = vector.min()
        if max_value is None:
            max_value = vector.max()
        _vec = vector[(vector >= min_value)*(vector <= max_value)]
        _vec.sort()
        binning = [min_value]
        for i in range(1, num_bins):
            index = int(i*len(_vec)/float(num_bins))
            binning.append(_vec[index])
        binning.append(max_value)
        return numpy.array(binning)

    def process_kwargs(self):
        """Check the keyword arguments.
        """
        if self.get('outfile') is None:
            suffx = self.__class__.__name__.replace('xEventBinning', '').lower()
            evfile = self.event_file.file_path()
            outfile = evfile.replace('.fits', '_%s.fits' % suffx)
            self.set('outfile', outfile)

    def bin_(self):
        """Do-nothing method to be reimplemented in the derived classes.
        """
        pass


class xBinnedFileBase:

    """Base class for binned files.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        assert(file_path.endswith('.fits'))
        logger.info('Opening input binned file %s...' % file_path)
        self.hdu_list = fits.open(file_path)
        self.hdu_list.info()

    def primary_header(self):
        """Return the primary header.
        """
        return self.hdu_list[0].header

    def primary_header_keyword(self, key):
        """Return the value of a primary header keyword.
        """
        return self.hdu_list[0].header[key]

    def plot(self, *args, **kwargs):
        """Do nothing method.
        """
        logger.info('%s.plot() not implemented, yet.' % self.__class__.__name__)


class xBinTableHDUPHA1(xBinTableHDUBase):

    """Binary table for binned PHA1 data.
    """

    NAME = 'SPECTRUM'
    HEADER_KEYWORDS = [
        ('HDUCLASS', 'OGIP'),
        ('HDUCLAS1', 'SPECTRUM'),
        ('HDUCLAS2', 'TOTAL'),
        ('HDUCLAS3', 'RATE'),
        ('CHANTYPE', 'PI'),
        ('HDUVERS' , '1.2.1', 'OGIP version number'),
        ('TLMIN1'  , 0      , 'first channel number'),
        ('TLMAX1'  , 4095    , 'last channel number'),
        ('CORRSCAL', 1.     , 'scaling for correction file'),
        ('POISSERR', False  , 'use statistical errors'),
        ('BACKFILE', None),
        ('CORRFILE', None),
        ('SYS_ERR' , 0.),
        ('AREASCAL', 1.),
        ('BACKSCAL', 1.)
    ]
    DATA_SPECS = [
        ('CHANNEL' , 'J'),
        ('RATE'    , 'E', 'counts/s'),
        ('STAT_ERR', 'E', 'counts/s'),
    ]


class xEventBinningPHA1(xEventBinningBase):

    """Class for PHA1 binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        num_chans = evt_header['DETCHANS']
        total_time = self.event_file.total_good_time()
        binning = numpy.linspace(-0.5, num_chans -0.5, num_chans+1)
        n, bins = numpy.histogram(self.event_data['PHA'], bins=binning)
        primary_hdu = self.build_primary_hdu()
        data = [numpy.arange(num_chans),
                n/total_time,
                numpy.sqrt(n)/total_time
        ]
        spec_hdu = xBinTableHDUPHA1(data)
        spec_hdu.setup_header(self.event_file.primary_keywords())
        irf_name = evt_header['IRFNAME']
        keywords = [('EXPOSURE', total_time, 'exposure time'),
                    ('RESPFILE', irf_file_path(irf_name, 'rmf')),
                    ('ANCRFILE', irf_file_path(irf_name, 'arf'))]
        spec_hdu.setup_header(keywords)
        hdu_list = fits.HDUList([primary_hdu, spec_hdu])
        hdu_list.info()
        logger.info('Writing binned PHA1 data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedCountSpectrum(xBinnedFileBase):

    """Binned phasogram.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xBinnedFileBase.__init__(self, file_path)
        self.data = self.hdu_list['SPECTRUM'].data
        self.channel = self.data['CHANNEL']
        self.rate = self.data['RATE']
        self.error = self.data['STAT_ERR']

    def plot(self, show=True):
        """Overloaded plot method.
        """
        plt.errorbar(self.channel, self.rate, yerr=self.error, fmt='o')
        plt.xlabel('PHA')
        plt.ylabel('Rate [Hz]')
        plt.yscale('log')
        if show:
            plt.show()


class xEventBinningCMAP(xEventBinningBase):

    """Class for CMAP binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        primary_header = self.event_file.hdu_list['PRIMARY'].header
        if self.get('xref') is None:
            self.set('xref', primary_header['ROIRA'])
        if self.get('yref') is None:
            self.set('yref', primary_header['ROIDEC'])

    def bin_(self):
        """Overloaded method.
        """
        if self.get('mc'):
            ra = self.event_data['MC_RA']
            dec = self.event_data['MC_DEC']
        else:
            ra = self.event_data['RA']
            dec = self.event_data['DEC']
        xref = self.get('xref')
        yref = self.get('yref')
        nxpix = self.get('nxpix')
        nypix = self.get('nypix')
        pixsize = self.get('binsz')/3600.
        proj = self.get('proj')
        sidex = nxpix*pixsize
        sidey = nypix*pixsize
        logger.info('Output image dimensions are %.1f x %.1f arcmin.' %\
                    (sidex*60, sidey*60))
        binsx = numpy.linspace(0, nxpix, nxpix + 1)
        binsy = numpy.linspace(0, nypix, nypix + 1)
        # Build the WCS object
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0.5*nxpix, 0.5*nypix]
        w.wcs.cdelt = [-pixsize, pixsize]
        w.wcs.crval = [xref, yref]
        w.wcs.ctype = ['RA---%s' % proj, 'DEC--%s' % proj]
        w.wcs.equinox = 2000.
        w.wcs.radesys = 'ICRS'
        header = w.to_header()
        # And here we need to tweak the header by hand to replicate what we
        # do in xEventBinningBase.build_primary_hdu() for the other binning
        # algorithms.
        header.set('BINALG', self.get('algorithm'),'the binning algorithm used')
        for key, val, comment in self.event_file.primary_keywords():
            header.set(key, val, comment)
        header['COMMENT'] = '%s run with kwargs %s' %\
                            (self.__class__.__name__, self.kwargs)
        # Ready to go!
        pix = w.wcs_world2pix(zip(ra, dec), 1)
        n, x, y = numpy.histogram2d(pix[:,1], pix[:,0], bins=(binsx, binsy))
        hdu = fits.PrimaryHDU(n, header=header)
        logger.info('Writing binned CMAP data to %s...' % self.get('outfile'))
        hdu.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedMap:

    """Display interface to binned CMAP files.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        self.image = xFITSImage(file_path, build_cdf=False)

    def plot(self, show=True,subplot=(1,1,1)):
        """Plot the data.
        """
        return self.image.plot(show=show,subplot=subplot)


class xEventBinningARATE(xEventBinningBase):

    """Class for ARATE binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        primary_header = self.event_file.hdu_list['PRIMARY'].header
        if self.get('xref') is None:
            self.set('xref', primary_header['ROIRA'])
        if self.get('yref') is None:
            self.set('yref', primary_header['ROIDEC'])

    def bin_(self):
        """Overloaded method.

        Warning
        -------
        This is horrible, as I shamelessly copied the CMAP code and added the
        missing bits. We should refactor the code in common instead.
        """
        if self.get('mc'):
            ra = self.event_data['MC_RA']
            dec = self.event_data['MC_DEC']
        else:
            ra = self.event_data['RA']
            dec = self.event_data['DEC']
        xref = self.get('xref')
        yref = self.get('yref')
        nxpix = self.get('nxpix')
        nypix = self.get('nypix')
        pixsize = self.get('binsz')/3600.
        proj = self.get('proj')
        sidex = nxpix*pixsize
        sidey = nypix*pixsize
        logger.info('Output image dimensions are %.1f x %.1f arcmin.' %\
                    (sidex*60, sidey*60))
        binsx = numpy.linspace(0, nxpix, nxpix + 1)
        binsy = numpy.linspace(0, nypix, nypix + 1)
        # Build the WCS object
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0.5*nxpix, 0.5*nypix]
        w.wcs.cdelt = [-pixsize, pixsize]
        w.wcs.crval = [xref, yref]
        w.wcs.ctype = ['RA---%s' % proj, 'DEC--%s' % proj]
        w.wcs.equinox = 2000.
        w.wcs.radesys = 'ICRS'
        header = w.to_header()
        # And here we need to tweak the header by hand to replicate what we
        # do in xEventBinningBase.build_primary_hdu() for the other binning
        # algorithms.
        header.set('BINALG', self.get('algorithm'),'the binning algorithm used')
        for key, val, comment in self.event_file.primary_keywords():
            header.set(key, val, comment)
        header['COMMENT'] = '%s run with kwargs %s' %\
                            (self.__class__.__name__, self.kwargs)
        # Ready to go!
        pix = w.wcs_world2pix(zip(ra, dec), 1)
        n, x, y = numpy.histogram2d(pix[:,1], pix[:,0], bins=(binsx, binsy))
        # Divide by the time to get the rate in counts per pixel.
        n /= self.event_file.total_good_time()
        # Ok, self.get('binsz') is the pixel size in arcseconds, and we have
        # to convert a solid angle in the sky to an area onto the detector.
        # self.get('focalscale') is the linear scale from arcmin in the sky to
        # mm on the focal plane.
        n /= (self.get('focalscale')*self.get('binsz')/60.)**2
        # And now, normalize to the number of telescopes on the focal plane.
        n /= self.get('ntelescopes')
        hdu = fits.PrimaryHDU(n, header=header)
        logger.info('Writing binned CMAP data to %s...' % self.get('outfile'))
        hdu.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedAreaRateMap:

    """Display interface to binned ARATE files.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        self.image = xFITSImage(file_path, build_cdf=False)

    def plot(self, show=True, subplot=(1,1,1)):
        """Plot the data.
        """
        return self.image.plot(show=show, zlabel='Counts / mm2 / s / GPD',
                               subplot=subplot)

    

class xBinTableHDULC(xBinTableHDUBase):

    """Binary table for binned LC data.
    """

    NAME = 'RATE'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('TIME'   , 'D', 's'     , 'time of the bin center'),
        ('TIMEDEL', 'D', 's'     , 'bin size'),
        ('COUNTS' , 'J', 'counts', 'photon counts'),
        ('ERROR'  , 'E', 'counts', 'statistical errors')
    ]


class xEventBinningLC(xEventBinningBase):

    """Class for LC binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        if self.get('tstart') is None:
            self.set('tstart', self.event_file.min_good_time())
        if self.get('tstop') is None:
            self.set('tstop', self.event_file.max_good_time())

    def make_binning(self):
        """Build the light-curve binning.
        """
        tbinalg = self.get('tbinalg')
        tstart = self.get('tstart')
        tstop = self.get('tstop')
        tbins = self.get('tbins')
        if tbinalg == 'LIN':
            return numpy.linspace(tstart, tstop, tbins + 1)
        elif tbinalg == 'LOG':
            return numpy.logspace(numpy.log10(tstart), numpy.log10(tstop),
                                  tbins + 1)
        elif tbinalg == 'FILE':
            tbinfile = self.get('tbinfile')
            assert tbinfile is not None
            return self.read_binning(tbinfile)
        abort('tbinalg %s not implemented yet' % tbinalg)

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        counts, edges = numpy.histogram(self.event_data['TIME'],
                                        bins=self.make_binning())
        primary_hdu = self.build_primary_hdu()
        data = [self.bin_centers(edges),
                self.bin_widths(edges),
                counts,
                numpy.sqrt(counts)
        ]
        rate_hdu = xBinTableHDULC(data)
        rate_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, rate_hdu, gti_hdu])
        hdu_list.info()
        logger.info('Writing binned LC data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedLightCurve(xBinnedFileBase):

    """Binned light curve.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xBinnedFileBase.__init__(self, file_path)
        self.data = self.hdu_list['RATE'].data
        self.time = self.data['TIME']
        self.timedel = self.data['TIMEDEL']
        self.counts = self.data['COUNTS']
        self.error = self.data['ERROR']

    def plot(self, show=True):
        """Overloaded plot method.
        """
        plt.errorbar(self.time, self.counts/self.timedel,
                     yerr=self.error/self.timedel, fmt='o')
        plt.xlabel('Time [s]')
        plt.ylabel('Rate [Hz]')
        if show:
            plt.show()


class xBinTableHDUPHASG(xBinTableHDUBase):

    """Binary table for binned PHASG data.
    """

    NAME = 'RATE'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('PHASE'   , 'D', 's'     , 'phase of the bin center'),
        ('PHASEDEL', 'D', 's'     , 'phase bin size'),
        ('COUNTS'  , 'J', 'counts', 'photon counts'),
        ('ERROR'   , 'E', 'counts', 'statistical errors')
    ]


class xEventBinningPHASG(xEventBinningBase):

    """Class for LC binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)

    def make_binning(self):
        """Build the light-curve binning.
        """
        phasebins = self.get('phasebins')
        return numpy.linspace(0., 1., phasebins + 1)

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        counts, edges = numpy.histogram(self.event_data['PHASE'],
                                        bins=self.make_binning())
        primary_hdu = self.build_primary_hdu()
        data = [self.bin_centers(edges),
                self.bin_widths(edges),
                counts,
                numpy.sqrt(counts)
        ]
        rate_hdu = xBinTableHDUPHASG(data)
        rate_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, rate_hdu, gti_hdu])
        hdu_list.info()
        logger.info('Writing binned PHASG data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedPhasogram(xBinnedFileBase):

    """Binned phasogram.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xBinnedFileBase.__init__(self, file_path)
        self.data = self.hdu_list['RATE'].data
        self.phase = self.data['PHASE']
        self.phase_delta = self.data['PHASEDEL']
        self.counts = self.data['COUNTS']
        self.error = self.data['ERROR']

    def plot(self, show=True, **kwargs):
        """Overloaded plot method.
        """
        if not kwargs.has_key('fmt'):
            kwargs['fmt'] = 'o'
        plt.errorbar(self.phase, self.counts, yerr=self.error,
                     xerr=0.5*self.phase_delta, **kwargs)
        plt.xlabel('Phase')
        plt.ylabel('Counts/bin')
        if show:
            plt.show()


class xBinTableHDUMCUBE(xBinTableHDUBase):

    """Binary table for binned MCUBE data.

    Mind the field for the actual phi distribution depends on the binning,
    which is specified at run time, and therefore the corresponding data specs
    must be set dinamically.
    """

    NAME = 'MODULATION'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('ENERGY_LO'   , 'E', 'keV'),
        ('ENERGY_HI'   , 'E', 'keV'),
        ('ENERGY_MEAN' , 'E', 'keV'),
        ('EFFECTIVE_MU', 'E'),
        ('COUNTS'      , 'J'),
        ('MDP 99%'     , 'E')
    ]

    @classmethod
    def set_phi_spec(cls, phibins):
        """Add the specification for the PHIHIST field.
        """
        phi_specs = ('PHI_HIST', '%dJ' % phibins)
        if 'PHI_HIST' in [spec[0] for spec in cls.DATA_SPECS]:
            cls.DATA_SPECS[-1] = phi_specs
        else:
            cls.DATA_SPECS.append(phi_specs)


class xEventBinningMCUBE(xEventBinningBase):

    """Class for MCUBE binning.
    """

    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)

    def make_binning(self):
        """Build the modulation cube binning.
        """
        ebinalg = self.get('ebinalg')
        emin = self.get('emin')
        emax = self.get('emax')
        ebins = self.get('ebins')
        if ebinalg == 'LIN':
            ebinning = numpy.linspace(emin, emax, ebins + 1)
        elif ebinalg == 'LOG':
            ebinning = numpy.linspace(numpy.log10(emin), numpy.log10(emax),
                                      ebins + 1)
        elif ebinalg == 'EQP':
            if self.get('mc'):
                energy = self.event_data['MC_ENERGY']
            else:
                energy = self.event_data['ENERGY']
            ebinning = self.equipopulated_binning(ebins, energy, emin, emax)
        elif ebinalg == 'FILE':
            ebinfile = self.get('ebinfile')
            assert ebinfile is not None
            ebinning = self.read_binning(ebinfile)
        elif ebinalg == 'LIST':
            ebinning = self.get('ebinning')
            assert isinstance(ebinning, list)
            ebinning = numpy.array(ebinning, 'd')
        else:
            abort('ebinalg %s not implemented yet' % ebinalg)
        phibinning = numpy.linspace(0, 2*numpy.pi, self.get('phibins') + 1)
        return (ebinning, phibinning)

    def bin_(self):
        """Overloaded method.
        """
        from ximpol.irf import load_mrf
        from ximpol.irf.mrf import mdp99
        modf = load_mrf(self.event_file.irf_name())
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        if self.get('mc'):
            energy = self.event_data['MC_ENERGY']
        else:
            energy = self.event_data['ENERGY']
        phi = self.event_data['PE_ANGLE']
        phi_hist, xedges, yedges = numpy.histogram2d(energy, phi,
                                                     bins=self.make_binning())
        primary_hdu = self.build_primary_hdu()
        emin, emax = xedges[:-1], xedges[1:]
        emean = []
        effmu = []
        ncounts = []
        mdp = []
        # If there's more than one bin in energy, we're also interested in all
        # the basic quantities in the entire energy range. This involves
        # extending the emin and emax arrays and summing up all the
        # phi histograms across the various energy slices.
        if len(emin) > 1:
            emin = numpy.append(emin, emin[0])
            emax = numpy.append(emax, emax[-1])
            _d1, _d2 = phi_hist.shape
            phi_hist_sum = numpy.sum(phi_hist, axis=0).reshape((1, _d2))
            phi_hist = numpy.append(phi_hist, phi_hist_sum, axis=0)
        for _emin, _emax in zip(emin, emax):
            _mask = (energy > _emin)*(energy < _emax)
            _energy = energy[_mask]
            _emean = numpy.mean(_energy)
            _effmu = modf.weighted_average(_energy)
            _ncounts = len(_energy)
            _mdp = mdp99(_effmu, _ncounts)
            emean.append(_emean)
            effmu.append(_effmu)
            ncounts.append(_ncounts)
            mdp.append(_mdp)
        data = [emin, emax, emean, effmu, ncounts, mdp, phi_hist]
        xBinTableHDUMCUBE.set_phi_spec(self.get('phibins'))
        mcube_hdu = xBinTableHDUMCUBE(data)
        mcube_hdu.setup_header(self.event_file.primary_keywords())
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, mcube_hdu, gti_hdu])
        hdu_list.info()
        logger.info('Writing binned MCUBE data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinnedModulationCube(xBinnedFileBase):

    """Read-mode interface to a MCUBE FITS file.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xBinnedFileBase.__init__(self, file_path)
        self.data = self.hdu_list['MODULATION'].data
        self.emin = self.data['ENERGY_LO']
        self.emax = self.data['ENERGY_HI']
        self.emean = self.data['ENERGY_MEAN']
        self.phi_y = self.data['PHI_HIST']
        self.effective_mu = self.data['EFFECTIVE_MU']
        self.counts = self.data['COUNTS']
        self.mdp99 = self.data['MDP 99%']
        phibins = self.phi_y.shape[1]
        self.phi_binning = numpy.linspace(0, 2*numpy.pi, phibins + 1)
        self.phi_x = 0.5*(self.phi_binning[:-1] + self.phi_binning[1:])
        from ximpol.irf import load_mrf
        irf_name = self.hdu_list['PRIMARY'].header['IRFNAME']
        self.fit_results = []

    def __len__(self):
        """Return the number of energy bins.
        """
        return len(self.emin)

    def fit_bin(self, i):
        """Fit the azimuthal distribution for the i-th energy slice.
        """
        hist = (self.phi_y[i], self.phi_binning, None)
        _fit_results = xAzimuthalResponseGenerator.fit_histogram(hist)
        _fit_results.set_polarization(self.effective_mu[i])
        logger.info(_fit_results)
        self.fit_results.append(_fit_results)
        return _fit_results

    def fit(self):
        """Fit the azimuthal distribution for all the energy bins.
        """
        self.fit_results = []
        for i in range(len(self)):
            self.fit_bin(i)

    def plot_bin(self, i, show=True, fit=True):
        """Plot the azimuthal distribution for the i-th energy slice.
        """
        _emin = self.emin[i]
        _emax = self.emax[i]
        _emean = self.emean[i]
        label = '%.2f-%.2f $<$%.2f$>$ keV' % (_emin, _emax, _emean)
        plt.errorbar(self.phi_x, self.phi_y[i], yerr=numpy.sqrt(self.phi_y[i]),
                     fmt='o')
        if fit:
            fit_results = self.fit_bin(i)
            fit_results.plot(label=label)

        plt.axis([0., 2*numpy.pi, 0.0, 1.2*self.phi_y[i].max()])
        plt.xlabel('Azimuthal angle [rad]')
        plt.ylabel('Counts/bin')
        plt.text(0.02, 0.92, label, transform=plt.gca().transAxes,
                 fontsize=15)
        if show:
            plt.show()

    def plot_polarization_degree(self, show=True, **kwargs):
        """Plot the polarization degree as a function of energy.
        """
        if self.fit_results == []:
            self.fit()
        _x = self.emean
        _dx = [self.emean - self.emin, self.emax - self.emean]
        _y = [r.polarization_degree for r in self.fit_results]
        _dy = [r.polarization_degree_error for r in self.fit_results]
        # If there's more than one energy binning we also fit the entire
        # energy interval, but we don't want the corresponding data point to
        # appear in the plot, se we brutally get rid of it.
        if len(self.fit_results) > 1:
            _x = _x[:-1]
            _dx = [_x - self.emin[:-1], self.emax[:-1] - _x]
            _y = _y[:-1]
            _dy = _dy[:-1]
        fig = plt.errorbar(_x, _y, _dy, _dx, fmt='o', **kwargs)
        plt.xlabel('Energy [keV]')
        plt.ylabel('Polarization degree')
        if show:
            plt.show()
        return fig

    def plot_polarization_angle(self, show=True, degree=False, **kwargs):
        """Plot the polarization angle as a function of energy.
        """
        if self.fit_results == []:
            self.fit()
        _x = self.emean
        _dx = [self.emean - self.emin, self.emax - self.emean]
        if degree:
            _y = [numpy.degrees(r.phase) for r in self.fit_results]
            _dy = [numpy.degrees(r.phase_error) for r in self.fit_results]
        else:
            _y = [(r.phase) for r in self.fit_results]
            _dy = [(r.phase_error) for r in self.fit_results]
        # If there's more than one energy binning we also fit the entire
        # energy interval, but we don't want the corresponding data point to
        # appear in the plot, se we brutally get rid of it.
        if len(self.fit_results) > 1:
            _x = _x[:-1]
            _dx = [_x - self.emin[:-1], self.emax[:-1] - _x]
            _y = _y[:-1]
            _dy = _dy[:-1]
        fig = plt.errorbar(_x, _y, _dy, _dx, fmt='o', **kwargs)
        plt.xlabel('Energy [keV]')
        if degree:
            plt.ylabel('Polarization angle [$^\circ$]')
        else:
            plt.ylabel('Polarization angle [rad]')
        if show:
            plt.show()
        return fig

    def plot(self, show=True, fit=True, analyze=True, xsubplot=0,
             full_range=True, simple_stat=False):
        """Plot the azimuthal distributions for all the energy bins.
        """
        if analyze:
            fit = True
        if fit:
            self.fit_results = []
        for i, _emean in enumerate(self.emean):
            if  full_range is False and i == (len(self.emean) - 1):
                continue
            if xsubplot == 0:
                plt.figure()
            else:
                plt.subplot(len(self.emean) + int(full_range) - 1, xsubplot + 1,
                            (xsubplot + 1)*(i + 1))
            self.plot_bin(i, False, False)
            if fit:
                _res = self.fit_bin(i)
                _res.plot(color=xpColor(i), simple_stat=simple_stat)
        if analyze:
            fig = plt.figure('Polarization degree vs. energy')
            self.plot_polarization_degree(show=False)
            fig = plt.figure('Polarization angle vs. energy')
            self.plot_polarization_angle(show=False)
        if show:
            plt.show()

class xBinTableHDUSCUBE(xBinTableHDUBase):

    """Binary table for binned SCUBE data.
    """

    NAME = 'EBOUNDS'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('ENERGY_LO'   , 'E', 'keV'),
        ('ENERGY_HI'   , 'E', 'keV'),
        ('EFFECTIVE_MU', 'E')
    ]


class xEventBinningSCUBE(xEventBinningBase):
    
    """Class for SCUBE binning.
    """
    
    def process_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.process_kwargs(self)
        primary_header = self.event_file.hdu_list['PRIMARY'].header
        if self.get('xref') is None:
            self.set('xref', primary_header['ROIRA'])
        if self.get('yref') is None:
            self.set('yref', primary_header['ROIDEC'])
    
    def make_binning(self):
        """Build the stokes cube energy binning.
        
        This is copied form modulation cube binning and should be factored out.
        """
        ebinalg = self.get('ebinalg')
        emin = self.get('emin')
        emax = self.get('emax')
        ebins = self.get('ebins')
        if ebinalg == 'LIN':
            ebinning = numpy.linspace(emin, emax, ebins + 1)
        elif ebinalg == 'LOG':
            ebinning = numpy.linspace(numpy.log10(emin), numpy.log10(emax),
                                      ebins + 1)
        elif ebinalg == 'EQP':
            if self.get('mc'):
                energy = self.event_data['MC_ENERGY']
            else:
                energy = self.event_data['ENERGY']
            ebinning = self.equipopulated_binning(ebins, energy, emin, emax)
        elif ebinalg == 'FILE':
            ebinfile = self.get('ebinfile')
            assert ebinfile is not None
            ebinning = self.read_binning(ebinfile)
        elif ebinalg == 'LIST':
            ebinning = self.get('ebinning')
            assert isinstance(ebinning, list)
            ebinning = numpy.array(ebinning, 'd')
        else:
            abort('ebinalg %s not implemented yet' % ebinalg)
        return ebinning
        
    def bin_(self):
        """Overloaded method.
        """
        from ximpol.irf import load_mrf
        modf = load_mrf(self.event_file.irf_name())
        if self.get('mc'):
            ra = self.event_data['MC_RA']
            dec = self.event_data['MC_DEC']
            energy = self.event_data['MC_ENERGY']
        else:
            ra = self.event_data['RA']
            dec = self.event_data['DEC']
            energy = self.event_data['ENERGY']
        phi = self.event_data['PE_ANGLE']
        # Calculate the stokes parameters
        I = numpy.ones(len(phi))
        Q = numpy.cos(2*phi)
        U = numpy.sin(2*phi)
        xref = self.get('xref')
        yref = self.get('yref')
        nxpix = self.get('nxpix')
        nypix = self.get('nypix')
        pixsize = self.get('binsz')/3600.
        proj = self.get('proj')
        sidex = nxpix*pixsize
        sidey = nypix*pixsize
        logger.info('Output image dimensions are %.1f x %.1f arcmin.' %\
                    (sidex*60, sidey*60))
        binsx = numpy.linspace(0, nxpix, nxpix + 1)
        binsy = numpy.linspace(0, nypix, nypix + 1)
        # Build the WCS object
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [0.5*nxpix, 0.5*nypix]
        w.wcs.cdelt = [-pixsize, pixsize]
        w.wcs.crval = [xref, yref]
        w.wcs.ctype = ['RA---%s' % proj, 'DEC--%s' % proj]
        w.wcs.equinox = 2000.
        w.wcs.radesys = 'ICRS'
        header = w.to_header()
        # And here we need to tweak the header by hand to replicate what we
        # do in xEventBinningBase.build_primary_hdu() for the other binning
        # algorithms.
        header.set('BINALG', self.get('algorithm'),'the binning algorithm used')
        for key, val, comment in self.event_file.primary_keywords():
            header.set(key, val, comment)
        header['COMMENT'] = '%s run with kwargs %s' %\
                            (self.__class__.__name__, self.kwargs)
        hdul = fits.HDUList()
        ebinning = self.make_binning()
        logger.info('Energy binning: %s' %ebinning)
        emin, emax = ebinning[:-1], ebinning[1:]
        mu_list = []
        #if len(emin) > 1:
        #    emin = numpy.append(emin, emin[0])
        #    emax = numpy.append(emax, emax[-1])
        for _emin, _emax in zip(emin, emax):
            _emask = (energy > _emin)*(energy < _emax)
            array_list = [I[_emask], Q[_emask], U[_emask], modf(energy)[_emask]]
            header.set('EMIN', _emin, 'minimum bin energy')
            header.set('EMAX', _emax, 'maximum bin energy')
            # Ready to go!
            pix = w.wcs_world2pix(zip(ra[_emask], dec[_emask]), 1)
            levels = []
            logger.info('Creating SCUBE for %.2f--%.2f keV...' %(_emin, _emax))
            for array in array_list:
                n, x, y = numpy.histogram2d(pix[:,1], pix[:,0],
                                            bins=(binsx, binsy), weights=array)
                levels.append(n)
            eff_mu = levels[3].sum()/levels[0].sum()
            mu_list.append(eff_mu)
            _mask = (levels[0] > 0)
            levels[3][_mask] /= levels[0][_mask]
            levels[2][_mask] = 2*levels[2][_mask]/levels[3][_mask]
            levels[1][_mask] = 2*levels[1][_mask]/levels[3][_mask]
            image = fits.ImageHDU(numpy.array(levels[:3]), header=header,
                                  name='STOKES: %.1f--%.1f keV'%(_emin, _emax))
            
            hdul.append(image)
        logger.info('Creating EBOUNDS...')
        data = [emin, emax, numpy.array(mu_list)]
        ebounds = xBinTableHDUSCUBE(data)
        ebounds.setup_header(self.event_file.primary_keywords())
        hdul.append(ebounds)
        logger.info('Writing SCUBE data to %s...' % self.get('outfile'))
        hdul.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')

class xBinnedStokesCube(xBinnedFileBase):

    """Read-mode interface to a SCUBE FITS file.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        xBinnedFileBase.__init__(self, file_path)
        self.hdu_list = fits.open(file_path)
        self.wcs = wcs.WCS(self.hdu_list[0].header)
        self.data = [self.hdu_list[i].data for i in range(0,
                     len(self.hdu_list)-1)]
        self.emin = self.hdu_list['EBOUNDS'].data['ENERGY_LO']
        self.emax = self.hdu_list['EBOUNDS'].data['ENERGY_HI']
        self.eff_mu = self.hdu_list['EBOUNDS'].data['EFFECTIVE_MU']
        
    def I(self, ebin=0):
        """Return the I counts map (slice=0) in the selected ebin.
        """
        map = self.data[ebin]
        return map[0,:,:]
    
    def Q(self, ebin=0):
        """Return the Q map (slice=1) in the selected ebin.
        """
        map = self.data[ebin]
        return map[1,:,:]
        
    def U(self, ebin=0):
        """Return the U map (slice=2) in the selected ebin.
        """
        map = self.data[ebin]
        return map[2,:,:]
    
    def getIQU(self, ebin=0):
        """Return the I, Q, U maps in the selected ebin.
        """
        return self.I(ebin), self.Q(ebin), self.U(ebin)
    
    def smooth(self, ebin=0, width=1):
        """Smooth the I, Q, U maps making the sum over the adjacent
        pixels in the selected ebin.
        
        Warning
        -------
        The N outer pixels per side are set to zero (N==width).
        """
        I, Q, U = self.getIQU(ebin)
        I_new = numpy.zeros(I.shape)
        Q_new = numpy.zeros(Q.shape)
        U_new = numpy.zeros(U.shape)
        for i in range(width, len(I[:,0])-width):
            for j in range(width, len(I[0,:])-width):
                I_new[i,j] = I[i-width:i+width+1,j-width:j+width+1].sum()
                Q_new[i,j] = Q[i-width:i+width+1,j-width:j+width+1].sum()
                U_new[i,j] = U[i-width:i+width+1,j-width:j+width+1].sum()
        return I_new, Q_new, U_new
    
    def polarization_degree_angle(self, ebin, smooth=0, I_cut=15, sigma=2,
        degree=True):
        """Return the polarization degree, angle and errors maps in the
        selected ebin.
        """
        if smooth is 0:
            I, Q, U = self.getIQU(ebin)
        else:
            I, Q, U = self.smooth(ebin, width=smooth)
        eff_mu = self.eff_mu[ebin]
        _mask = (I >= I_cut)
        pol_deg = numpy.zeros(I.shape)
        pol_ang = numpy.zeros(I.shape)
        pol_deg_err = numpy.zeros(I.shape)
        pol_ang_err = numpy.zeros(I.shape)
        pol_deg[_mask] = numpy.sqrt(Q[_mask]**2+U[_mask]**2)/I[_mask]
        pol_deg_err[_mask] =\
            numpy.sqrt((2./eff_mu**2-pol_deg[_mask]**2)/(I[_mask]-1))
        pol_ang[_mask] = 0.5*numpy.arctan2(U[_mask], Q[_mask])
        pol_ang[pol_ang<0.] += numpy.pi
        pol_ang_err[_mask] =\
            1./(pol_deg[_mask]*eff_mu*numpy.sqrt(2*(I[_mask]-1)))
        if degree:
            pol_ang = numpy.degrees(pol_ang)
            pol_ang_err = numpy.degrees(pol_ang_err)
        _sigma_mask = (((pol_deg + sigma*pol_deg_err) > 1) +\
                       ((pol_deg - sigma*pol_deg_err) < 0))
        pol_deg[_sigma_mask] = 0.
        pol_deg_err[_sigma_mask] = 0.
        pol_ang[_sigma_mask] = 0.
        pol_ang_err[_sigma_mask] = 0.
        return [pol_deg, pol_deg_err, pol_ang, pol_ang_err]
    
    def plot_polarization_degree_angle(self, pol_list, show=True):
        """Plot the polarization degree and angle maps.
        """
        zlabel = 'Polarization degree'
        fig_deg, fig_deg_err = self._make_pol_plot(pol_list[0], pol_list[1],
                                                   zlabel, show=show)
        zlabel = 'Polarization angle'
        fig_ang, fig_ang_err = self._make_pol_plot(pol_list[2], pol_list[3],
                                                   zlabel, show=show)
        return [fig_deg, fig_deg_err, fig_ang, fig_ang_err]

    def _make_pol_plot(self, pol, pol_err, zlabel, show=True):
        """Make the polarization plot.
        """
        hdul = fits.HDUList()
        img_pol = fits.ImageHDU(pol, header=self.wcs.to_header())
        img_pol_err =fits.ImageHDU(pol_err, header=self.wcs.to_header())
        hdul.append(img_pol)
        hdul.append(img_pol_err)
        fig_pol = self._make_plot(hdul[0], zlabel, show=show)
        zlabel += ' error'
        fig_pol_err = self._make_plot(hdul[1], zlabel, show=show)
        return fig_pol, fig_pol_err
    
    def plot(self, ebin=0, slice=0, show=True, zlabel=None, subplot=(1, 1, 1)):
        """Plot the Stokes parameters (I, Q or U) map in the selected ebin.
        """
        label_list = ['I', 'Q', 'U']
        if zlabel is None:
            zlabel = label_list[slice]
        return self._make_plot(self.hdu_list[ebin], zlabel, slice=slice,
                               show=show, subplot=subplot)
    
    def _make_plot(self, hdul, zlabel, slice=0, show=True, subplot=(1, 1, 1)):
        """Make and show the plot.

        This is using aplpy to render the image.

        Warning
        -------
        We have to figure out the subplot issue, here. I put in a horrible
        hack to recover the previous behavior when there's only one
        subplot.
        """
        import aplpy
        from ximpol.utils.matplotlib_ import context_no_grids
        with context_no_grids():
            if subplot == (1, 1, 1):
                fig = aplpy.FITSFigure(hdul,figure=plt.figure(), slices=[slice])
            else:
                fig = aplpy.FITSFigure(hdul, slices=[slice], subplot=subplot,
                figure=plt.figure(0,figsize=(10*subplot[1], 10*subplot[0])))
        fig.add_grid()
        fig.show_colorscale(cmap = 'afmhot')#cmap='plasma')
        fig.add_colorbar()
        fig.colorbar.set_axis_label_text(zlabel)
        if show:
            plt.show()
        return fig
