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


from ximpol.irf.base import xColDefsBase


MATRIX_HEADER_SPECS = [
    ('TUNIT1'  , 'keV     ', 'physical unit of field'),
    ('TUNIT2'  , 'keV     ', 'physical unit of field'),
    ('EXTNAME' , 'MATRIX'  , 'name of this binary table extension'),
    ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
    ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
    ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
    ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
    ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the form'),
    ('DETCHANS', '1024'    , 'total number of detector channels'),
    ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
    ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility'),
    ('CHANTYPE' , 'PI '    , 'Detector Channel Type in use (PHA or PI)'),
    ('CREATOR' , None      , 's/w task which wrote this dataset'),
    ('RESPFILE', None      , None)
]


EBOUNDS_HEADER_SPECS = [
    ('TUNIT2'  , 'keV ', 'physical unit of field'),
    ('TUNIT3'  , 'keV ',  'physical unit of field'),
    ('EXTNAME' , 'EBOUNDS ',' Name of this binary table extension'),
    ('CHANTYPE', 'PI  ',' Channel type'),
    ('CONTENT' , 'Response Matrix' ,' File content'),
    ('HDUCLASS', 'OGIP    ' ,' Format conforms to OGIP/GSFC conventions'),
    ('HDUCLAS1', 'RESPONSE', 'Extension contains response data  '),
    ('HDUCLAS2', 'EBOUNDS ',' Extension contains EBOUNDS'),
    ('DETCHANS', 1024     , 'Total number of detector channels')
    ]



class xColDefsMATRIX(xColDefsBase):

    """
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E'),
        ('ENERG_HI', 'E'),
        ('N_GRP'   , 'I'),
        ('F_CHAN'  , 'I'),
        ('N_CHAN'  , 'I'),
        ('MATRIX'  , '1024E')
    ]


class xColDefsEBOUNDS(xColDefsBase):

    """
    """

    COLUMN_SPECS = [
        ('CHANNEL', 'I'),
        ('E_MIN'  , 'E'),
        ('E_MAX'  , 'E')
    ]
