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


"""TODO: clean up the header and get rid of the hard-coded stuff.
"""

SPECTRUM_HEADER_SPECS = [
    ('EXTNAME' , 'SPECTRUM', 'name of this binary table extension'),
    ('HDUCLASS', 'OGIP'    , None),
    ('HDUCLAS1', 'SPECTRUM', None),
    ('HDUCLAS2', 'TOTAL'   , None),
    ('HDUCLAS3', 'RATE'    , None),
    ('CHANTYPE', 'PI'      , None),
    ('HDUVERS' , '1.2.1'   , 'OGIP version number'),
    ('TELESCOP', None      , None),
    ('INSTRUME', None      , None),
    ('DETNAM'  , None      , None),
    ('FILTER'  , None      , None),
    ('DATAMODE', None      , None),
    ('DETCHANS', 256       , 'Number of channels in spectrum'),
    ('TLMIN1'  , 0         , 'First channel number'),
    ('EXPOSURE', 1.        , 'Exposure time'),
    ('CORRSCAL', 1.        , 'Scaling for correction file'),
    ('POISSERR', 'T'       , 'Is error Poisson ?'),
    ('RESPFILE', None      , None),
    ('ANCRFILE', None      , None),
    ('BACKFILE', None      , None),
    ('CORRFILE', None      , None),
    ('SYS_ERR' , 0.        , None),
    ('QUALITY' , 0         , None),
    ('GROUPING', 0         , None),
    ('AREASCAL', 1.        , None),
    ('BACKSCAL', 1.        , None)
]


class xColDefsSpectrum(xColDefsBase):

    """
    """

    COLUMN_SPECS = [
        ('CHANNEL' , 'J', None),
        ('RATE'    , 'E', 'counts/s'),
        ('STAT_ERR', 'E', None),
    ]
