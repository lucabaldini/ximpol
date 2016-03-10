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


"""Facilities related to FITS I/O.
"""


import time
import numpy

from  astropy.io import fits

from ximpol.__version__ import TAG


"""
From https://heasarc.gsfc.nasa.gov/docs/software/fitsio/c/c_user/node20.html
Codes for the data type of binary table columns and/or for the
data type of variables when reading or writing keywords or data:

                              DATATYPE               TFORM CODE
  #define TBIT          1  /*                            'X' */
  #define TBYTE        11  /* 8-bit unsigned byte,       'B' */
  #define TLOGICAL     14  /* logicals (int for keywords     */
                           /*  and char for table cols   'L' */
  #define TSTRING      16  /* ASCII string,              'A' */
  #define TSHORT       21  /* signed short,              'I' */
  #define TLONG        41  /* signed long,                   */
  #define TLONGLONG    81  /* 64-bit long signed integer 'K' */
  #define TFLOAT       42  /* single precision float,    'E' */
  #define TDOUBLE      82  /* double precision float,    'D' */
  #define TCOMPLEX     83  /* complex (pair of floats)   'C' */
  #define TDBLCOMPLEX 163  /* double complex (2 doubles) 'M' */

  The following data type codes are also supported by CFITSIO:
  #define TINT         31  /* int                            */
  #define TSBYTE       12  /* 8-bit signed byte,         'S' */
  #define TUINT        30  /* unsigned int               'V' */
  #define TUSHORT      20  /* unsigned short             'U'  */
  #define TULONG       40  /* unsigned long                  */

  The following data type code is only for use with fits\_get\_coltype
  #define TINT32BIT    41  /* signed 32-bit int,         'J' */
"""


FITS_TO_NUMPY_TYPE_DICT = {
    'E': numpy.float32,
    'D': numpy.float64,
    'I': numpy.int16,
    'J': numpy.int32
    }


class xHDUBase:

    """Base class for FITS HDU.
    """

    def add_comment(self, comment):
        """Add a comment to the table header.
        """
        self.header['COMMENT'] = comment

    def add_keyword(self, key, value, comment=''):
        """Add a keyword to the table header.
        """
        self.header.set(key, value, comment)

    def set_keyword_comment(self, keyword, comment):
        """Set the comment for a header keyword.
        """
        self.header.comments[keyword] = comment

    def setup_header(self, keywords=[], comments=[]):
        """Update the table header with arbitrary additional information.
        """
        for item in keywords:
            if len(item) == 3:
                key, value, comment = item
            elif len(item) == 2:
                key, value = item
                comment = ''
            self.add_keyword(key, value, comment)
        for comment in comments:
            self.add_comment(comment)

    def __str__(self):
        """String formatting.
        """
        return repr(self.header)


class xPrimaryHDU(fits.PrimaryHDU, xHDUBase):

    """Class describing a primary HDU to be written in a FITS file.

    This is initializing a standard astropy.io.fits.PrimaryHDU object and
    adding the creator and creation time fields.

    Arguments
    ---------
    creator : str
       The application that created the header (defaults to `ximpol`, and
       the ximpol tag is automatically added.)
    """

    def __init__(self, creator='ximpol', keywords=[], comments=[]):
        """Constructor.
        """
        fits.PrimaryHDU.__init__(self)
        self.header.set('CREATOR',
                        '%s %s' % (creator, TAG),
                        's/w task which wrote this dataset')
        self.header.set('DATE',
                        time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime()),
                        'file creation date (YYYY-MM-DDThh:mm:ss UT')
        self.setup_header(keywords, comments)


class xBinTableHDUBase(fits.BinTableHDU, xHDUBase):

    """Binary table HDU class.

    This is a small wrapper around a standard binary table to facilitate
    customizations.
    """

    NAME = None
    HEADER_KEYWORDS = []
    HEADER_COMMENTS = []
    DATA_SPECS = []

    def __init__(self, data=None, keywords=[], comments=[]):
        """
        """
        if data is not None:
            assert(len(data) == len(self.DATA_SPECS))
        cols = []
        _kwcomments = {}
        for i, item in enumerate(self.DATA_SPECS):
            if len(item) == 4:
                name, format_, units, comment = item
                _kwcomments[name] = comment
            elif len(item) == 3:
                name, format_, units = item
                description = ''
            elif len(item) == 2:
                name, format_ = item
                units = None
                description = ''
            if data is not None:
                col = fits.Column(name, format_, units, array=data[i])
            else:
                col = fits.Column(name, format_, units)
            cols.append(col)
        data = fits.FITS_rec.from_columns(cols)
        fits.BinTableHDU.__init__(self, data)
        # Set the extension name, if necessary.
        if self.NAME is not None:
            self.set_ext_name(self.NAME)
        # Add the additional keywords and comments to the header.
        self.setup_header(self.HEADER_KEYWORDS + keywords,
                          self.HEADER_COMMENTS + comments)
        # And we need to loop one more time to take care of the comments on the
        # columns, if any (could not find a way to do this on the columns
        # directly).
        for i, col in enumerate(self.columns):
            if _kwcomments.has_key(col.name):
                comment = _kwcomments[col.name]
                self.set_keyword_comment('TTYPE%d' % (i + 1), comment)

    @classmethod
    def spec_names(self):
        """Return the name of the data fields specified in the SPEC class
        member.
        """
        return [item[0] for item in self.DATA_SPECS]

    @classmethod
    def spec_names_and_types(self):
        """Return the name of the data fields specified in the SPEC class
        member.
        """
        return [item[0:2] for item in self.DATA_SPECS]

    def set_ext_name(self, name):
        """Set the extension name for the binary table.
        """
        self.add_keyword('EXTNAME', name, 'name of this binary table extension')

    def __str__(self):
        """String formatting.
        """
        return '%s\n%s' % (repr(self.header), self.data)



def main():
    """
    """
    import numpy

    class CustomBinTable(xBinTableHDUBase):

        NAME = 'CUSTOM'
        HEADER_KEYWORDS = [
            ('PKEY', None)
        ]
        HEADER_COMMENTS = [
            'A comment'
        ]
        DATA_SPECS = [
            ('ENERGY' , 'E', 'keV'  , 'The energy value'),
            ('CHANNEL', 'I'),
            ('AEFF'   , 'E', 'cm**2')
        ]

    hdu = xPrimaryHDU()
    print(hdu)
    n = 20
    data = [numpy.linspace(1., 10., n), numpy.arange(n), numpy.ones(n)]
    table1 = CustomBinTable()
    print(table1)
    keywords = [
        ('AKEY1', 0, 'A keyword'),
        ('AKEY2', 1)
    ]
    comments = ['Howdy, partner?']
    table2 = CustomBinTable(data, keywords, comments)
    print(table2)


if __name__ == '__main__':
    main()
