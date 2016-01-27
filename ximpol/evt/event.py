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
    SPECS = [
        ('TIME'    , 'D', 's'      , 'Event time in seconds'),
        ('PHA'     , 'I', None     , 'Uncorrected event channel'),
        ('PE_ANGLE', 'E', 'degrees', 'Reconstructed photoelectron angle'),
        ('ENERGY'  , 'E', 'KeV'    , 'Reconstructed event energy'),
        ('RA'      , 'E', 'degrees', 'Reconstructed right ascension'),
        ('DEC'     , 'E', 'degrees', 'Reconstructed declination')
    ]


class xBinTableHDUMonteCarloEvents(xBinTableHDUBase):

    """Binary table description for the EVENTS extension of the observation
    output files, including the additional Monte Carlo fields.
    """

    NAME = xBinTableHDUEvents.NAME
    SPECS = xBinTableHDUEvents.SPECS + [
        ('MC_ENERGY', 'E', 'keV'    , 'Monte Carlo event energy'),
        ('MC_RA'    , 'E', 'degrees', 'Monte Carlo right ascension'),
        ('MC_DEC'   , 'E', 'degrees', 'Monte Carlo declination'),
        ('MC_SRC_ID', 'I', None     , 'Monte Carlo souce identifier')
    ]


class xBinTableHDUGTI(xBinTableHDUBase):

    """Binary tablefor the good time intervals (GTI).
    """

    NAME = 'GTI'
    SPECS = [
        ('START', 'D', 's', 'GTI start time'),
        ('STOP' , 'D', 's', 'GTI stop time')
    ]


class xBinTableHDURoiTable(xBinTableHDUBase):

    """Binary tablefor the good time intervals (GTI).
    """

    NAME = 'ROITABLE'
    SPECS = [
        ('SRCID'  , 'I'  , None, 'Source identifier'),
        ('SRCNAME', 'A20', None, 'Source name')
    ]


class xMonteCarloEventList(dict):

    """Class describing a Monte Carlo event list.
    """

    def __init__(self):
        """Constructor.
        """
        for name in xBinTableHDUMonteCarloEvents.spec_names():
            self[name] = numpy.array([])
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
        if isinstance(data, numbers.Number):
            _col = numpy.empty(self.length)
            _col.fill(data)
            data = _col
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

    def write_fits(self, file_path, gti_list=[], roi_model=None):
        """Write the event list and associated ancillary information to file.

        Arguments
        ---------
        file_path : str
            The path to the output file.

        gti : list of 2-elements (start, stop) tuples
            The list of good time intervals.

        roi_model :
           The ROI model used to generate the event list.
        """
        primary_hdu = xPrimaryHDU()
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


def main():
    """
    """
    pass


if __name__ == '__main__':
    main()
