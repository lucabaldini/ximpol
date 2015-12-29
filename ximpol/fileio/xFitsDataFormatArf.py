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



from ximpol.fileio.xFitsDataFormatBase import xFitsDataFormatBase
from ximpol.__logging__ import logger



class xFitsDataFormatArf(xFitsDataFormatBase):

    """ Specification for the arf data format.
    """

    PRIMARY_HEADER_SPECS = [
        ('SIMPLE'  , True      , 'file does conform to FITS standard'),
        ('BITPIX'  , -32       , 'number of bits per data pixel'),
        ('NAXIS'   , 0         , 'number of data axes'),
        ('EXTEND'  , True      , 'FITS dataset may contain extensions'),
        ('DATE'    , None      , 'file creation date (YYYY-MM-DDThh:mm:ss UT'),
        ('CREATOR' , None      , 's/w task which wrote this dataset')
    ]

    SPECRESP_HEADER_SPECS = [
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
        ('TTYPE3'  , 'SPECRESP', 'label for field   3'),
        ('TFORM3'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT3'  , 'cm**2   ', 'physical unit of field'),
        ('EXTNAME' , 'SPECRESP', 'name of this binary table extension'),
        ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
        ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
        ('HDUCLAS2', 'SPECRESP', 'dataset contains spectral response'),
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

    SPECRESP_DATA_SPECS = [
        ('ENERG_LO', 'E'),
        ('ENERG_HI', 'E'),
        ('SPECRESP', 'E')
    ]
    
    @staticmethod
    def primaryHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatArf.PRIMARY_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)

    @staticmethod
    def specrespHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatArf.SPECRESP_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)

    @staticmethod
    def specrespColumns(arrays, **kwargs):
        """
        """
        specs = xFitsDataFormatArf.SPECRESP_DATA_SPECS
        return xFitsDataFormatBase.columns(specs, arrays, **kwargs)




def test():
    """ Test code.
    """
    logger.info('Creating .arf PRIMARY header...')
    p = xFitsDataFormatArf.primaryHeader()
    print(repr(p))
    logger.info('Creating .arf SPECRESP header...')
    s = xFitsDataFormatArf.specrespHeader(TELESCOP = 'XIPE', INSTRUME = 'GPD')
    print(repr(s))



if __name__ == '__main__':
    test()
    
    
