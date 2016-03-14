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


"""Module encapsulating the event structure and related facilities.
"""


import numpy
import numbers
from astropy.io import fits

from ximpol.utils.logging_ import logger
from ximpol.core.fitsio import xPrimaryHDU, xBinTableHDUBase
from ximpol.core.fitsio import FITS_TO_NUMPY_TYPE_DICT


class xBinTableHDUEvents(xBinTableHDUBase):

    """Binary table description for the EVENTS extension of the observation
    output files.

    This is partially modeled based on the `XMM-Newton data format
    <http://www.mpe.mpg.de/xray/wave/xmm/cookbook/EPIC_PN/event_descr.php>`_.

    A recap of the basic FITS data types is `here
    <https://pythonhosted.org/pyfits/users_guide/users_table.html>`_.

    This bit is separated from the additional Monte Carlo event fields, as
    in the future the "real" part of the data format might be stored
    elsewhere.
    """

    NAME = 'EVENTS'
    DATA_SPECS = [
        ('TIME'    , 'D', 's'      , 'event time in seconds'),
        ('PHA'     , 'I', None     , 'uncorrected event channel'),
        ('PE_ANGLE', 'E', 'degrees', 'reconstructed photoelectron angle'),
        ('ENERGY'  , 'E', 'KeV'    , 'reconstructed event energy'),
        ('RA'      , 'E', 'degrees', 'reconstructed right ascension'),
        ('DEC'     , 'E', 'degrees', 'reconstructed declination')
    ]


class xBinTableHDUMonteCarloEvents(xBinTableHDUBase):

    """Binary table description for the EVENTS extension of the observation
    output files, including the additional Monte Carlo fields.
    """

    NAME = xBinTableHDUEvents.NAME
    DATA_SPECS = xBinTableHDUEvents.DATA_SPECS + [
        ('PHASE'    , 'E', None     , 'Event phase (for periodic sources)'),
        ('MC_ENERGY', 'E', 'keV'    , 'Monte Carlo event energy'),
        ('MC_RA'    , 'E', 'degrees', 'Monte Carlo right ascension'),
        ('MC_DEC'   , 'E', 'degrees', 'Monte Carlo declination'),
        ('MC_SRC_ID', 'I', None     , 'Monte Carlo souce identifier')
    ]


class xBinTableHDUGTI(xBinTableHDUBase):

    """Binary tablefor the good time intervals (GTI).
    """

    NAME = 'GTI'
    DATA_SPECS = [
        ('START', 'D', 's', 'GTI start time'),
        ('STOP' , 'D', 's', 'GTI stop time')
    ]


class xBinTableHDURoiTable(xBinTableHDUBase):

    """Binary tablefor the good time intervals (GTI).
    """

    NAME = 'ROITABLE'
    DATA_SPECS = [
        ('SRCID'  , 'I'  , None, 'source identifier'),
        ('SRCNAME', 'A20', None, 'source name')
    ]


class xMonteCarloEventList(dict):

    """Class describing a Monte Carlo event list.
    """

    def __init__(self):
        """Constructor.
        """
        for name, dtype in xBinTableHDUMonteCarloEvents.spec_names_and_types():
            dtype = FITS_TO_NUMPY_TYPE_DICT[dtype]
            self[name] = numpy.array([], dtype)
        self.length = 0

    def __len__(self):
        """Return the length of the event list.
        """
        return self.length

    def set_column(self, name, data):
        """Set a column array.

        Arguments
        ---------
        name : string
            The name of the column

        data : array or number
            The actual data to put in the column. (If `data` is a number,\
            an array of the proper length is automatically created, assuming\
            that the `lenght` class member is defined.)
        """
        assert self.has_key(name)
        dtype = self[name].dtype
        if isinstance(data, numbers.Number):
            data = numpy.full(self.length, data, dtype)
        if self.length > 0:
            assert(len(data) == self.length)
        else:
            self.length = len(data)
        self[name] = data

    def __add__(self, other):
        """Concatenate two event lists.
        """
        _list = xMonteCarloEventList()
        for name in xBinTableHDUMonteCarloEvents.spec_names():
            _list.set_column(name,numpy.append(self[name], other[name]))
        return _list

    def sort(self):
        """Sort the event list based on the event time.
        """
        _index = numpy.argsort(self['TIME'])
        for name in xBinTableHDUMonteCarloEvents.spec_names():
            self.set_column(name, self[name][_index])

    def write_fits(self, file_path, simulation_info):
        """Write the event list and associated ancillary information to file.

        Arguments
        ---------
        file_path : str
            The path to the output file.

        simulation_info :
            A generic container with all the relevant information about the
            simulation.

        Warning
        -------
        The information about the detector and telescope should be in the
        primary header of the IRF tables, and that's where we should be
        retrieving it from. (See issue #49.)
        """
        primary_hdu = xPrimaryHDU()
        roi_model = simulation_info.roi_model
        irf_name = simulation_info.irf_name
        ebounds_header = simulation_info.edisp.hdu_list['EBOUNDS'].header
        gti_list = simulation_info.gti_list
        keywords = [
            ('ROIRA'   , roi_model.ra , 'right ascension of the ROI center'),
            ('ROIDEC'  , roi_model.dec, 'declination of the ROI center'),
            ('EQUINOX' , 2000.        , 'equinox for RA and DEC'),
            ('IRFNAME' , irf_name     , 'name of the IRFs used for the MC'),
            ('TELESCOP', ebounds_header['TELESCOP']),
            ('INSTRUME', ebounds_header['INSTRUME']),
            ('DETNAM'  , ebounds_header['DETNAM']),
            ('DETCHANS', ebounds_header['DETCHANS'])
        ]
        primary_hdu.setup_header(keywords)
        data = [self[name] for name in\
                xBinTableHDUMonteCarloEvents.spec_names()]
        event_hdu = xBinTableHDUMonteCarloEvents(data)
        _start = numpy.array([gti[0] for gti in gti_list])
        _stop = numpy.array([gti[1] for gti in gti_list])
        gti_hdu = xBinTableHDUGTI([_start, _stop])
        _src_id = numpy.array([src.identifier for src in roi_model.values()])
        _src_name = numpy.array([src.name for src in roi_model.values()])
        roi_hdu = xBinTableHDURoiTable([_src_id, _src_name])
        hdu_list = fits.HDUList([primary_hdu, event_hdu, gti_hdu, roi_hdu])
        hdu_list.info()
        hdu_list.writeto(file_path, clobber=True)
        logger.info('Event list written to %s...' % file_path)


class xEventFile:

    """Read-mode interface to event files.
    """

    def __init__(self, file_path):
        """Constructor.
        """
        assert(file_path.endswith('.fits'))
        logger.info('Opening input event file %s...' % file_path)
        self.hdu_list = fits.open(file_path)
        self.hdu_list.info()
        self.event_data = self.hdu_list['EVENTS'].data
        self.roi_table = self.build_roi_table()

    def close(self):
        """Close the HDU list.
        """
        self.hdu_list.close()

    def num_events(self):
        """Return the total number of events in the event file.
        """
        return len(self.event_data)

    def file_path(self):
        """Return the path to the underlying file.
        """
        return self.hdu_list.filename()

    def primary_keywords(self):
        """Return a list of all the relevant keywords in the
        PRIMARY header.

        This is used, e.g.,  to propagate all the necessary information (such as
        the ROI and the IRFs used) from the event files to the binned
        data files that are created from them.
        """
        _header = self.hdu_list['PRIMARY'].header
        keywords = []
        for key in ['ROIRA', 'ROIDEC', 'EQUINOX', 'IRFNAME', 'TELESCOP',
                    'INSTRUME', 'DETCHANS']:
            keywords.append((key, _header[key], _header.comments[key]))
        return keywords

    def build_roi_table(self):
        """Rebuild the ROI table based in the information in the ROITABLE
        extension of the event file.
        """
        roi_table = {}
        _src_id = self.hdu_list['ROITABLE'].data['SRCID']
        _src_name = self.hdu_list['ROITABLE'].data['SRCNAME']
        for _id, _name in zip(_src_id, _src_name):
            roi_table[_name] = _id
        return roi_table

    def roi_center(self):
        """Return the righ ascension and declination of the center of the
        ROI, as written in the header of the primary extension.
        """
        return self.hdu_list['PRIMARY'].header['ROIRA'],\
            self.hdu_list['PRIMARY'].header['ROIDEC']

    def irf_name(self):
        """Return the name of the IRF set used to run the simulation.
        """
        return self.hdu_list['PRIMARY'].header['IRFNAME']

    def total_good_time(self):
        """Return the sum of all the GTIs in the event file.
        """
        _start = self.hdu_list['GTI'].data['START']
        _stop = self.hdu_list['GTI'].data['STOP']
        return (_stop - _start).sum()

    def min_good_time(self):
        """Return the smaller START time for the GTIs in the event file.
        """
        return self.hdu_list['GTI'].data['START'].min()

    def max_good_time(self):
        """Return the largest STOP time for the GTIs in the event file.
        """
        return self.hdu_list['GTI'].data['STOP'].max()

    def source_id(self, source_name):
        """Return the source identifier for a given source name in the ROI.
        """
        return self.roi_table[source_name]
