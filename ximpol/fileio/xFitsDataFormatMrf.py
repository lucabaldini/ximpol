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


from ximpol.fileio.xFitsDataFormatBase import xFitsDataFormatBase
from ximpol.__logging__ import logger


class xFitsDataFormatMrf(xFitsDataFormatBase):

    """Specification for the arf data format.
    """

    PRIMARY_HEADER_SPECS = [
        ('SIMPLE'  , True      , 'file does conform to FITS standard'),
        ('BITPIX'  , -32       , 'number of bits per data pixel'),
        ('NAXIS'   , 0         , 'number of data axes'),
        ('EXTEND'  , True      , 'FITS dataset may contain extensions'),
        ('DATE'    , None      , 'file creation date (YYYY-MM-DDThh:mm:ss UT'),
        ('CREATOR' , None      , 's/w task which wrote this dataset')
    ]

    MODFRESP_HEADER_SPECS = [
        ('XTENSION', 'BINTABLE', 'binary table extension'),
        ('BITPIX'  , 8         , '8-bit bytes'),
        ('NAXIS'   , 2         , '2-dimensional binary table'),
        ('NAXIS1'  , 12        , 'width of table in bytes'),
        ('NAXIS2'  , None      , 'number of rows in table'),
        ('PCOUNT'  , 0         , 'size of special data area'),
        ('GCOUNT'  , 1         , 'one data group (required keyword)'),
        ('TFIELDS' , 3         , 'number of fields in each row'),
        ('TTYPE1'  , 'ENERG_LO', 'label for field   1'),
        ('TFORM1'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT1'  , 'keV     ', 'physical unit of field'),
        ('TTYPE2'  , 'ENERG_HI', 'label for field   2'),
        ('TFORM2'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT2'  , 'keV     ', 'physical unit of field'),
        ('TTYPE3'  , 'MODFRESP', 'label for field   3'),
        ('TFORM3'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT3'  , None      , 'physical unit of field'),
        ('EXTNAME' , 'MODFRESP', 'name of this binary table extension'),
        ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
        ('HDUCLAS1', 'RESPONSE', 'dataset relates to polarimetric response'),
        ('HDUCLAS2', 'MODFRESP', 'dataset contains polarimetric response'),
        ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
        ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the forma'),
        ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
        ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility'),
        ('TELESCOP', None      , 'mission/satellite name'),
        ('INSTRUME', None      , 'instrument/detector name'),
        ('DETNAM'  , None      , 'specific detector name in use'),
        ('FILTER'  , None      , 'filter in use'),
        ('ARFVERSN', '1992a   ', 'Obsolete - included for backwards compatibility'),
        ('PHAFILE' , 'UNKNOWN ', 'PHA file for which this ARF created'),
        ('CREATOR' , None      , 's/w task which wrote this dataset'),
        ('RESPFILE', None      , None)
    ]

    MODFRESP_DATA_SPECS = [
        ('ENERG_LO', 'E'),
        ('ENERG_HI', 'E'),
        ('MODFRESP', 'E')
    ]

    @staticmethod
    def primaryHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatMrf.PRIMARY_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)

    @staticmethod
    def modfrespHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatMrf.MODFRESP_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)

    @staticmethod
    def modfrespColumns(arrays, **kwargs):
        """
        """
        specs = xFitsDataFormatMrf.MODFRESP_DATA_SPECS
        return xFitsDataFormatBase.columns(specs, arrays, **kwargs)


def main():
    logger.info('Creating .mrf PRIMARY header...')
    p = xFitsDataFormatMrf.primaryHeader()
    print(repr(p))
    logger.info('Creating .mrf MODFRESP header...')
    s = xFitsDataFormatMrf.modfrespHeader(TELESCOP = 'XIPE', INSTRUME = 'GPD')
    print(repr(s))


if __name__ == '__main__':
    main()
