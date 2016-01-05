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


"""Basic format definitions related to the instrument response functions.
"""


import time

from  astropy.io import fits

from ximpol.__version__ import TAG


"""Specifications for the FITS headers related to the OGIP standards.
(These are in common to all the fits files.)
"""

OGIP_HEADER_SPECS = [
    ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
    ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
    ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the forma'),
    ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
    ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility')
]



class xPrimaryHDU(fits.PrimaryHDU):

    """Class describing a primary HDU to be written in a FITS file.

    This is initializing a standard astropy.io.fits.PrimaryHDU object and
    adding the creator and creation time fields. (The constructor takes no
    arguments.)

    Examples
    --------
    >>> primary_hdu = xPrimaryHDU()
    >>> print(repr(primary_hdu.header))
    """

    def __init__(self):
        """Constructor.
        """
        fits.PrimaryHDU.__init__(self)
        creator = 'ximpol %s' % TAG
        date = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
        self.header.set('CREATOR', creator,
                        's/w task which wrote this dataset')
        self.header.set('DATE', date,
                        'file creation date (YYYY-MM-DDThh:mm:ss UT')



class xColDefsBase(fits.ColDefs):

    """Column definitions base class.

    This is a do-nothing class. Subclasses must define the COLUMN_SPECS
    class member in the form of a list of three-elements tuples containing
    (in this order):

    * the name of the column;
    * the format of the column content;
    * the physical units of the column, e.g.:

    >>> COLUMN_SPECS = [
    >>>       ('ENERG_LO', 'E', 'keV'),
    >>>       ('ENERG_HI', 'E', 'keV'),
    >>>       ('SPECRESP', 'E', 'cm**2')
    >>> ]

    The actual data are filled in at creation time.

    Args
    ----
    data : list of arrays
        The actual data to be put in the columns. (The number and shapes of the
        arrays have to match the column definitions.)
    """

    COLUMN_SPECS = []

    def __init__(self, data):
        """
        """
        cols = []
        for i, (name, fmt, units) in enumerate(self.COLUMN_SPECS):
            col = fits.Column(name, fmt, units, array = data[i])
            cols.append(col)
        fits.ColDefs.__init__(self, cols)
