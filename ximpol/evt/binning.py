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
from ximpol.utils.matplotlib_ import pyplot as plt

"""TODO: clean up the header and get rid of the hard-coded stuff.

We want to (loosely) model this on
http://fermi.gsfc.nasa.gov/ssc/data/analysis/scitools/help/gtbin.txt
"""


class xEventBinningBase:

    """Base class for the event binning.
    """

    VALID_KWARGS = [
        'outfile'
    ]

    def __init__(self, file_path, **kwargs):
        """Constructor.
        """
        self.event_file = xEventFile(file_path)
        self.event_data = self.event_file.event_data
        self.kwargs = kwargs
        self.check_kwargs()

    def check_kwargs(self):
        """Check the keyword arguments.
        """
        for key in self.kwargs.keys():
            if not key in self.VALID_KWARGS:
                abort('Invalid keyword argument "%s" for %s' %\
                      (key, self.__class__.__name__))
        if not 'outfile' in self.kwargs.keys():
            _sfx = self.__class__.__name__.replace('xEventBinning', '').lower()
            _in_path = self.event_file.file_path()
            _out_path = _in_path.replace('.fits', '_%s.fits' % _sfx)
            logger.info('outfile not specified, using %s...' % _out_path)
            self.kwargs['outfile'] = _out_path

    def bin_(self):
        """Do the actual binning.
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

    VALID_KWARGS = xEventBinningBase.VALID_KWARGS

    def check_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.check_kwargs(self)

    def bin_(self):
        """Overloaded method.
        """
        evt_header = self.event_file.hdu_list['PRIMARY'].header
        num_chans = evt_header['DETCHANS']
        good_time = self.event_file.good_time()
        binning = numpy.linspace(-0.5, num_chans - 0.5, num_chans)
        n, bins, patches = plt.hist(self.event_data['PHA'], bins=binning)
        primary_hdu = xPrimaryHDU()
        data = [numpy.arange(num_chans),
                n/good_time,
                numpy.sqrt(n)/good_time
        ]
        spec_hdu = xBinTableHDUPHA1(data)
        keywords = [
            ('TELESCOP', evt_header['TELESCOP']),
            ('INSTRUME', evt_header['INSTRUME']),
            ('DETCHANS', num_chans, 'number of channels in spectrum'),
            ('EXPOSURE', good_time, 'exposure time'),
        ]
        spec_hdu.setup_header(keywords)
        hdu_list = fits.HDUList([primary_hdu, spec_hdu])
        hdu_list.info()
        logger.info('Writing binned (PHA1) data to %s...' % \
                    self.kwargs['outfile'])
        hdu_list.writeto(self.kwargs['outfile'], clobber=True)
        logger.info('Done.')


class xEventBinningCMAP(xEventBinningBase):

    """Class for CMAP binning.
    """

    VALID_KWARGS = xEventBinningBase.VALID_KWARGS

    def check_kwargs(self):
        """Overloaded method.
        """
        xEventBinningBase.check_kwargs(self)

    def bin_(self):
        """Overloaded method.

        Warning
        -------
        This needs to be completely revised.
        """
        mc = False
        nside = 256
        side = 5./60
        primary_header = self.event_file.hdu_list['PRIMARY'].header
        event_data = self.event_file.hdu_list['EVENTS'].data
        if mc:
            ra = event_data['MC_RA']
            dec = event_data['MC_DEC']
        else:
            ra = event_data['RA']
            dec = event_data['DEC']
        roi_ra = primary_header['ROIRA']
        roi_dec = primary_header['ROIDEC']
        delta = side/nside
        binning = numpy.linspace(0, nside, nside + 1)
        # Build the WCS object
        w = wcs.WCS(naxis=2)
        w.wcs.crpix = [nside, 0.]
        w.wcs.cdelt = [-delta, delta]
        w.wcs.crval = [roi_ra - 0.5*side, roi_dec - 0.5*side]
        w.wcs.ctype = ["RA---TAN", "DEC--TAN"]
        w.wcs.equinox = 2000
        header = w.to_header()
        pix = w.wcs_world2pix(zip(ra, dec), 1)
        n, x, y = numpy.histogram2d(pix[:,1], pix[:,0], bins=(binning, binning))
        hdu = fits.PrimaryHDU(n, header=header)
        hdu.writeto(self.kwargs['outfile'], clobber=True)
        logger.info('Writing binned (CMAP) data to %s...' %\
                    self.kwargs['outfile'])
