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


import time

from  astropy.io import fits

from ximpol.__version__ import TAG


class xPrimaryHDU(fits.PrimaryHDU):

    """
    """

    def __init__(self):
        """
        """
        fits.PrimaryHDU.__init__(self)
        creator = 'ximpol %s' % TAG
        date = time.strftime('%Y-%m-%dT%H:%M:%S', time.gmtime())
        self.header.set('CREATOR', creator,
                        's/w task which wrote this dataset')
        self.header.set('DATE', date,
                        'file creation date (YYYY-MM-DDThh:mm:ss UT')



class xColDefsBase(fits.ColDefs):

    """
    """

    COLUMN_SPECS = []

    def __init__(self, columns):
        """
        """
        cols = []
        for i, (name, fmt) in enumerate(self.COLUMN_SPECS):
            col = fits.Column(name = name, format = fmt, array = columns[i])
            cols.append(col)
        fits.ColDefs.__init__(self, cols)
