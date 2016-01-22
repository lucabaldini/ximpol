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
from ximpol.irf.base import xColDefsBase
from ximpol.irf.base import xPrimaryHDU, update_header


EVENT_HEADER_SPECS = [
    ('EXTNAME' , 'EVENTS', 'name of this binary table extension')
]


class xColDefsEvent(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for an event FITS files.

    This is partially modeled based on the `XMM-Newton data format
    <http://www.mpe.mpg.de/xray/wave/xmm/cookbook/EPIC_PN/event_descr.php>`_.

    A recap of the basic FITS data types is `here
    <https://pythonhosted.org/pyfits/users_guide/users_table.html>`_.

    This bit is separated from the additional Monte Carlo event fields, as
    in the future the "real" part of the data format might be stored
    elsewhere.
    """

    COLUMN_SPECS = [
        ('TIME'    , 'E', 's'),
        ('PHA'     , 'I', None),
        #('DETX'    , 'E', 'mm'),
        #('DETY'    , 'E', 'mm'),
        #('X'       , 'E', 'arcsecs'),
        #('Y'       , 'E', 'arcsecs'),
        ('PE_ANGLE', 'E', 'degrees'),
        # And, for convenience, we add these three, too.
        ('ENERGY'  , 'E', 'KeV'),
        ('RA'      , 'E', 'degrees'),
        ('DEC'     , 'E', 'degrees')
    ]


class xColDefsMonteCarloEvent(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for a Monte Carlo event FITS files.
    """

    COLUMN_SPECS = xColDefsEvent.COLUMN_SPECS + [
        ('MC_ENERGY', 'E', 'keV'),
        ('MC_RA'    , 'E', 'degrees'),
        ('MC_DEC'   , 'E', 'degrees'),
        ('MC_SRC_ID', 'I', None)
    ]


class xMonteCarloEventList(dict):

    """Class describing a Monte Carlo event list.
    """

    def __init__(self):
        """Constructor.
        """
        for (name, fmt, units) in xColDefsMonteCarloEvent.COLUMN_SPECS:
            self[name] = None
        self.length = None

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
        if self.length is not None:
            assert(len(data) == self.length)
        else:
            self.length = len(data)
        self[name] = data

    def __add__(self, other):
        """Concatenate two event lists.
        """
        _list = xMonteCarloEventList()
        for (name, fmt, units) in xColDefsMonteCarloEvent.COLUMN_SPECS:
            _list.set_column(numpy.concatenate(self[name], other[name]))
        return _list

    def sort(self):
        """Sort an event list.
        """
        _index = numpy.argsort(self['TIME'])
        for (name, fmt, units) in xColDefsMonteCarloEvent.COLUMN_SPECS:
            self.set_column(self[name][_index])

    def write_fits(self, file_path):
        """Write the event list to file.
        """
        primary_hdu = xPrimaryHDU()
        data = []
        for (name, fmt, units) in xColDefsMonteCarloEvent.COLUMN_SPECS:
            data.append(self[name])
        cols = xColDefsMonteCarloEvent(data)
        event_hdu = fits.BinTableHDU.from_columns(cols)
        update_header(event_hdu, EVENT_HEADER_SPECS)
        hdulist = fits.HDUList([primary_hdu, event_hdu])
        hdulist.info()
        hdulist.writeto(file_path, clobber=True)
        logger.info('Event list written to %s...' % file_path)


def main():
    """
    """
    pass


if __name__ == '__main__':
    main()
