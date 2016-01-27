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

from  astropy.io import fits

from ximpol.__version__ import TAG


class xPrimaryHDU(fits.PrimaryHDU):

    """Class describing a primary HDU to be written in a FITS file.

    This is initializing a standard astropy.io.fits.PrimaryHDU object and
    adding the creator and creation time fields.

    Arguments
    ---------
    creator : str
       The application that created the header (defaults to `ximpol`, and
       the ximpol tag is automatically added.)
    """

    def __init__(self, creator='ximpol'):
        """Constructor.
        """
        fits.PrimaryHDU.__init__(self)
        self.header.set('CREATOR',
                        '%s %s' % (creator, TAG),
                        's/w task which wrote this dataset')
        self.header.set('DATE',
                        time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime()),
                        'file creation date (YYYY-MM-DDThh:mm:ss UT')

    def __str__(self):
        """String formatting.
        """
        return repr(self.header)


class xBinTableHDUBase(fits.BinTableHDU):

    """Binary table HDU class.

    This is a small wrapper around a standard binary table to facilitate
    customizations.
    """

    NAME = None
    SPECS = []

    def __init__(self, data=None, keywords=[], comments=[]):
        """
        """
        if data is not None:
            assert(len(data) == len(self.SPECS))
        cols = []
        _kwcomments = {}
        for i, item in enumerate(self.SPECS):
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
        self.setup_header(keywords, comments)
        # And we need to loop one more time to take care of the comments on the
        # columns, if any (could not find a way to do this on the columns
        # directly).
        for i, col in enumerate(self.columns):
            if _kwcomments.has_key(col.name):
                self.set_header_comment('TTYPE%d' % (i + 1),
                                        _kwcomments[col.name])

    @classmethod
    def spec_names(self):
        """Return the name of the data fields specified in the SPEC class
        member.
        """
        return [item[0] for item in self.SPECS]

    def set_ext_name(self, name):
        """Set the extension name for the binary table.
        """
        self.add_header_keyword('EXTNAME', name,
                                'name of this binary table extension')

    def add_header_keyword(self, key, value, comment=''):
        """Add a keyword to the table header.
        """
        self.header.set(key, value, comment)

    def add_header_comment(self, comment):
        """Add a comment to the table header.
        """
        self.header['COMMENT'] = comment

    def set_header_comment(self, keyword, comment):
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
            self.add_header_keyword(key, value, comment)
        for comment in comments:
            self.add_header_comment(comment)

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
        SPECS = [
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
        ('KEY1', 0, 'A keyword'),
        ('KEY2', 1)
    ]
    comments = ['Howdy, partner?']
    table2 = CustomBinTable(data, keywords, comments)
    print(table2)


if __name__ == '__main__':
    main()
