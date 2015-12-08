#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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



import time
from astropy.io import fits

from ximpol.__version__ import TAG



class xFitsDataFormatBase:

    """ Basic interface to a generic set of specifications for writing out
    fits files.
    """

    @staticmethod
    def header(specs, comments = [], **kwargs):
        """ Create a fits header based on a given set of format specifications.
        """
        header = fits.Header()
        if not kwargs.has_key('CREATOR'):
            kwargs['CREATOR'] = 'ximpol %s' % TAG
        if not kwargs.has_key('DATE'):
            kwargs['DATE'] = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
        for key, default, comment in specs:
            try:
                value = kwargs[key]
            except KeyError:
                if default is None:
                    default = 'NONE'
                value = default
            if comment is not None:
                header[key] = (value, comment)
            else:
                header[key] = value
        for comment in comments:
            header['COMMENT'] = comment
        return header



def test():
    """ Test code.
    """
    specs = [
        ('SIMPLE'  , True      , 'file does conform to FITS standard'),
        ('BITPIX'  , -32       , 'number of bits per data pixel'),
        ('NAXIS'   , 0         , 'number of data axes'),
        ('EXTEND'  , True      , 'FITS dataset may contain extensions'),
        ('DATE'    , None      , 'file creation date (YYYY-MM-DDThh:mm:ss UT'),
        ('CREATOR' , None      , 's/w task which wrote this dataset')
    ]
    comments = [
        'This is a test file.',
        'And cool, too.'
    ]
    h = xFitsDataFormatBase.header(specs, comments, CREATOR = 'Luca')
    print(repr(h))



if __name__ == '__main__':
    test()
