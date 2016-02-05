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


"""Module encapsulating the FITS spectra structure and related facilities.
"""

import numpy
from astropy.io import fits
from astropy import wcs

from ximpol.utils.logging_ import logger, abort
from ximpol.evt.event import xEventFile
from ximpol.core.fitsio import xBinTableHDUBase
from ximpol.irf.base import xColDefsBase
from ximpol.irf.base import xPrimaryHDU, update_header


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

    @classmethod
    def bin_centers(self, bin_edges):
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
    def bin_widths(self, bin_edges):
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
        ('CORRSCAL', 1.     , 'scaling for correction file'),
        ('POISSERR', 'T'    , 'is error Poisson?'),
        ('RESPFILE', None),
        ('ANCRFILE', None),
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
        binning = numpy.linspace(-0.5, num_chans - 0.5, num_chans)
        n, bins = numpy.histogram(self.event_data['PHA'], bins=binning)
        primary_hdu = xPrimaryHDU()
        data = [numpy.arange(num_chans),
                n/total_time,
                numpy.sqrt(n)/total_time
        ]
        spec_hdu = xBinTableHDUPHA1(data)
        keywords = [
            ('TELESCOP', evt_header['TELESCOP']),
            ('INSTRUME', evt_header['INSTRUME']),
            ('DETCHANS', num_chans, 'number of channels in spectrum'),
            ('EXPOSURE', total_time, 'exposure time'),
        ]
        spec_hdu.setup_header(keywords)
        hdu_list = fits.HDUList([primary_hdu, spec_hdu])
        hdu_list.info()
        logger.info('Writing binned PHA1 data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


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
        w.wcs.crpix = [nxpix, 0.]
        w.wcs.cdelt = [-pixsize, pixsize]
        w.wcs.crval = [xref - 0.5*sidex, yref - 0.5*sidey]
        w.wcs.ctype = ['RA---%s' % proj, 'DEC--%s' % proj]
        w.wcs.equinox = 2000
        header = w.to_header()
        pix = w.wcs_world2pix(zip(ra, dec), 1)
        n, x, y = numpy.histogram2d(pix[:,1], pix[:,0], bins=(binsx, binsy))
        hdu = fits.PrimaryHDU(n, header=header)
        logger.info('Writing binned CMAP data to %s...' % self.get('outfile'))
        hdu.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


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
        if self.get('tbinalg') == 'LIN':
            return numpy.linspace(self.get('tstart'),
                                  self.get('tstop'),
                                  self.get('tbins') + 1)
        if self.get('tbinalg') == 'LOG':
            return numpy.linspace(numpy.log10(self.get('tstart')),
                                  numpy.log10(self.get('tstop')),
                                  self.get('tbins') + 1)
        abort('%s not implemented yet' % self.get('tbinalg'))

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        counts, edges = numpy.histogram(self.event_data['TIME'],
                                        bins=self.make_binning())
        primary_hdu = xPrimaryHDU()
        data = [self.bin_centers(edges),
                self.bin_widths(edges),
                counts,
                numpy.sqrt(counts)
        ]
        rate_hdu = xBinTableHDULC(data)
        keywords = [
            ('TELESCOP', evt_header['TELESCOP']),
            ('INSTRUME', evt_header['INSTRUME'])
        ]
        rate_hdu.setup_header(keywords)
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, rate_hdu, gti_hdu])
        hdu_list.info()
        logger.info('Writing binned LC data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')


class xBinTableHDUMCUBE(xBinTableHDUBase):

    """Binary table for binned MCUBE data.

    Mind the field for the actual phi distribution depends on the binning,
    which is specified at run time, and therefore the corresponding data specs
    must be set dinamically.
    """

    NAME = 'MODULATION'
    HEADER_KEYWORDS = []
    DATA_SPECS = [
        ('ENERGY_LO'  , 'E', 'keV'),
        ('ENERGY_HI'  , 'E', 'keV'),
        ('ENERGY_MEAN', 'E', 'keV')
    ]

    @classmethod
    def add_phi_spec(self, phibins):
        """Add the specification for the PHIHIST field.
        """
        phi_specs = ('PHIHIST', '%dI' % phibins)
        self.DATA_SPECS.append(phi_specs)


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
        if self.get('ebinalg') == 'LIN':
            ebinning = numpy.linspace(self.get('emin'),
                                      self.get('emax'),
                                      self.get('ebins') + 1)
        elif self.get('ebinalg') == 'LOG':
            ebinning = numpy.linspace(numpy.log10(self.get('emin')),
                                      numpy.log10(self.get('emax')),
                                      self.get('ebins') + 1)
        else:
            abort('%s not implemented yet' % self.get('ebinalg'))
        phibinning = numpy.linspace(0, 2*numpy.pi, self.get('phibins') + 1)
        return (ebinning, phibinning)

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        if self.get('mc'):
            energy = self.event_data['MC_ENERGY']
        else:
            energy = self.event_data['ENERGY']
        phi = self.event_data['PE_ANGLE']
        counts, xedges, yedges = numpy.histogram2d(energy, phi,
                                                   bins=self.make_binning())
        primary_hdu = xPrimaryHDU()
        emin, emax = xedges[:-1], xedges[1:]
        emean = []
        for _emin, _emax in zip(emin, emax):
            emean.append(numpy.mean(energy[(energy > _emin)*(energy < _emax)]))
        data = [emin, emax, emean, counts]
        xBinTableHDUMCUBE.add_phi_spec(self.get('phibins'))
        mcube_hdu = xBinTableHDUMCUBE(data)
        keywords = [
            ('TELESCOP', evt_header['TELESCOP']),
            ('INSTRUME', evt_header['INSTRUME'])
        ]
        mcube_hdu.setup_header(keywords)
        gti_hdu = self.event_file.hdu_list['GTI']
        hdu_list = fits.HDUList([primary_hdu, mcube_hdu, gti_hdu])
        hdu_list.info()
        logger.info('Writing binned MCUBE data to %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')
