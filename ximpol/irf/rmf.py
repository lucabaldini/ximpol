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


from ximpol.irf.base import xColDefsBase, OGIP_HEADER_SPECS


"""Header specifications for the MATRIX extension of .rmf FITS files.
"""
MATRIX_HEADER_SPECS = [
    ('EXTNAME' , 'MATRIX'    , 'name of this binary table extension'),
    ('HDUCLAS1', 'RESPONSE'  , 'dataset relates to spectral response'),
    ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
    ('DETCHANS', '1024'      , 'total number of detector channels'),
    ('CHANTYPE', 'PI '       , 'Detector Channel Type in use (PHA or PI)')
] + OGIP_HEADER_SPECS


"""Header specifications for the EBOUNDS extension of .rmf FITS files.
"""
EBOUNDS_HEADER_SPECS = [
    ('EXTNAME' , 'EBOUNDS '  , 'Name of this binary table extension'),
    ('CHANTYPE', 'PI'        , 'Channel type'),
    ('CONTENT' , 'Response Matrix', 'File content'),
    ('HDUCLAS1', 'RESPONSE'  , 'Extension contains response data  '),
    ('HDUCLAS2', 'EBOUNDS '  , 'Extension contains EBOUNDS'),
    ('DETCHANS', 1024        , 'Total number of detector channels')
] + OGIP_HEADER_SPECS



class xColDefsMATRIX(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the MATRIX extension
    of .rmf FITS files.
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('N_GRP'   , 'I', None),
        ('F_CHAN'  , 'I', None),
        ('N_CHAN'  , 'I', None),
        ('MATRIX'  , '1024E', None)
    ]


class xColDefsEBOUNDS(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the EBOUNDS extension
    of .rmf FITS files.
    """

    COLUMN_SPECS = [
        ('CHANNEL', 'I', None),
        ('E_MIN'  , 'E', 'keV'),
        ('E_MAX'  , 'E', 'keV')
    ]
