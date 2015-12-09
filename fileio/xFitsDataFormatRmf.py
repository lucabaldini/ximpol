#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# * Copyright (C) 2015 Melissa Pesce-Rollins                          *
# *              (melissa.pesce.rollins@pi.infn.it)                   *
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



class xFitsDataFormatRmf(xFitsDataFormatBase):

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

    MODFRESP_HEADER_SPECS = [
        ('XTENSION', 'BINTABLE', 'binary table extension'),
        ('BITPIX'  , 8         , '8-bit bytes'),
        ('NAXIS'   , 2         , '2-dimensional binary table'),
        ('NAXIS1'  , 4110      , 'length of dimension 1'),
        ('NAXIS2'  , 2400      , 'length of dimension 2'),
        ('PCOUNT'  , 0         , 'umber of group parameters'),
        ('GCOUNT'  , 1         , 'number of groups'),
        ('TFIELDS' , 6         , 'number of table fields'),
        ('TTYPE1'  , 'ENERG_LO', 'label for field   1'),
        ('TFORM1'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT1'  , 'keV     ', 'physical unit of field'),
        ('TTYPE2'  , 'ENERG_HI', 'label for field   2'),
        ('TFORM2'  , 'E       ', 'data format of field: 4-byte REAL'),
        ('TUNIT2'  , 'keV     ', 'physical unit of field'),
        ('TTYPE3 ' , 'N_GRP '  , 'label for field 3: 2-byte INTEGER'),
        ('TFORM3 ' , 'I '      , 'total number of channel subsets'),
        ('TTYPE4 ' , 'F_CHAN ' , 'label for field 4:  2 or 4-byte INTEGER array'),
        ('TFORM4 ' , 'I'       , 'channel number of the start of each channel subset '),
        ('TTYPE5 ' , 'N_CHAN ' , 'label for field 5: 2 or 4-byte INTEGER vector'),
        ('TFORM5 ' , 'I '      , 'number of channels within each "channel subset" for the energy bin.'),
        ('TTYPE6 ' , 'MATRIX ' , 'label for field 6: REAL vector array, each element within is 4-byte'),
        ('TFORM6'  , '1024E '  , 'the response values for each channel subset for the energy bin.'),
        ('EXTNAME' , 'MATRIX'  , 'name of this binary table extension'),
        ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
        ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
        ('HDUCLAS2', 'RSP_MATRIX', 'dataset is a spectral response matrix'),
        ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
        ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the form'),
        ('DETCHANS', '1024'    , 'total number of detector channels'),
        ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
        ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility'),
        ('TELESCOP', None      , 'mission/satellite name'),
        ('INSTRUME', None      , 'instrument/detector name'),
        ('DETNAM'  , None      , 'specific detector name in use'),
        ('FILTER'  , None      , 'filter in use'),
        ('ARFVERSN', '1992a   ', 'Obsolete - included for backwards compatibility'),
        ('CHANTYPE' , 'PI '    , 'Detector Channel Type in use (PHA or PI)'),
        ('CREATOR' , None      , 's/w task which wrote this dataset'),
        ('RESPFILE', None      , None)
    ]

    MODFRESP_DATA_SPECS = [
        ('ENERG_LO', 'E'),
        ('ENERG_HI', 'E'),
        ('N_GRP'   , 'I'),
        ('F_CHAN'  , 'I'),
        ('N_CHAN'  , 'I'),
        ('MATRIX'  , '1024E')
    ]

    @staticmethod
    def primaryHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatRmf.PRIMARY_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)
    
    @staticmethod
    def modfrespHeader(comments = [], **kwargs):
        """
        """
        specs = xFitsDataFormatRmf.MODFRESP_HEADER_SPECS
        return xFitsDataFormatBase.header(specs, comments, **kwargs)

    @staticmethod
    def modfrespColumns(arrays, **kwargs):
        """
        """
        specs = xFitsDataFormatRmf.MODFRESP_DATA_SPECS
        return xFitsDataFormatBase.columns(specs, arrays, **kwargs)




def test():
    """ Test code.
    """
    logger.info('Creating .rmf PRIMARY header...')
    p = xFitsDataFormatRmf.primaryHeader()
    print(repr(p))
    logger.info('Creating .rmf MODFRESP header...')
    s = xFitsDataFormatRmf.modfrespHeader(TELESCOP = 'XIPE', INSTRUME = 'GPD')
    print(repr(s))



if __name__ == '__main__':
    test()
    
    
