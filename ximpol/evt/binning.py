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

from ximpol.utils.logging_ import logger
from ximpol.irf.base import xColDefsBase
from ximpol.irf.base import xPrimaryHDU, update_header


SPECTRUM_HEADER_SPECS = [
    ('EXTNAME' , 'SPECTRUM', 'name of this binary table extension')
]

class xColDefsSpectrum(xColDefsBase):

    """
    """

    COLUMN_SPECS = [
        ('CHANNEL' , 'J', None),
        ('RATE'    , 'E', 'counts/s'),
        ('STAT_ERR', 'E', None),
    ]
