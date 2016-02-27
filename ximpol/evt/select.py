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


import numpy
import time
from astropy.io import fits
from astropy.coordinates import SkyCoord
from astropy import wcs

from ximpol.utils.logging_ import logger, abort
from ximpol.evt.event import xEventFile
from ximpol.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ximpol.utils.matplotlib_ import pyplot as plt


class xEventSelect:

    """Base class for event subselection.
    """

    def __init__(self, file_path, **kwargs):
        """Constructor.
        """
        self.event_file = xEventFile(file_path)
        self.event_data = self.event_file.event_data
        self.kwargs = kwargs
        self.process_kwargs()

    def process_kwargs(self):
        """Check the keyword arguments.
        """
        if self.get('outfile') is None:
            evfile = self.event_file.file_path()
            outfile = evfile.replace('.fits', '_select.fits')
            self.set('outfile', outfile)
        if self.get('ra') is None:
            self.set('ra', self.event_file.roi_center()[0])
        if self.get('dec') is None:
            self.set('dec', self.event_file.roi_center()[1])

    def get(self, key, default=None):
        """Convenience method to address the keyword aguments.
        """
        return self.kwargs.get(key, default)

    def set(self, key, value):
        """Convenience method to set keyword arguments.
        """
        logger.info('Setting %s to %s...' % (key, value))
        self.kwargs[key] = value

    def select(self):
        """Select the events and write the output file.
        """
        logger.info('Running event selection with kwargs %s...' % self.kwargs)
        if self.get('mc'):
            evt_ra = self.event_data['MC_RA']
            evt_dec = self.event_data['MC_DEC']
            evt_energy = self.event_data['MC_ENERGY']
        else:
            evt_ra = self.event_data['RA']
            evt_dec = self.event_data['DEC']
            evt_energy = self.event_data['ENERGY']
        num_events = self.event_file.num_events()
        mask = numpy.ones(num_events, 'bool')
        if self.get('rad') is not None:
            evt_skycoord = SkyCoord(evt_ra, evt_dec, unit='deg')
            ref_skyccord = SkyCoord(self.get('ra'), self.get('dec'), unit='deg')
            separation = evt_skycoord.separation(ref_skyccord).arcmin
            mask *= (separation < self.get('rad'))
        if self.get('emin') is not None:
            mask *= (evt_energy > self.get('emin'))
        if self.get('emax') is not None:
            mask *= (evt_energy < self.get('emax'))
        if self.get('tmin') is not None:
            mask *= (self.event_data['TIME'] > self.get('tmin'))
        if self.get('tmax') is not None:
            mask *= (self.event_data['TIME'] < self.get('tmax'))
        if self.get('phasemin') is not None:
            mask *= (self.event_data['PHASE'] > self.get('phasemin'))
        if self.get('phasemax') is not None:
            mask *= (self.event_data['PHASE'] < self.get('phasemax'))
        for srcid in self.get('mcsrcid'):
            mask *= (self.event_data['MC_SRC_ID'] == srcid)
        events_hdu = self.event_file.hdu_list['EVENTS']
        events_hdu.data = events_hdu.data[mask]
        logger.info('Done, %d out of %d remaining...' %\
                    (len(events_hdu.data), num_events))
        primary_hdu = self.event_file.hdu_list['PRIMARY']
        timestamp = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
        comment = 'xEventSelect run on %s with kwargs %s' %\
                  (timestamp, self.kwargs)
        primary_hdu.header['COMMENT'] = comment
        gti_hdu = self.event_file.hdu_list['GTI']
        roi_hdu = self.event_file.hdu_list['ROITABLE']
        hdu_list = fits.HDUList([primary_hdu, events_hdu, gti_hdu, roi_hdu])
        hdu_list.info()
        logger.info('Writing data subselectionto %s...' % self.get('outfile'))
        hdu_list.writeto(self.get('outfile'), clobber=True)
        logger.info('Done.')
