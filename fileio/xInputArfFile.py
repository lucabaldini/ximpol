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



from astropy.io import fits

from ximpol.fileio.xInputFileBase import xInputFileBase
from ximpol.__logging__ import logger



class xInputArfFile(xInputFileBase):

    """
    """

    def __init__(self, filePath):
        """
        """
        xInputFileBase.__init__(self, filePath, 'arf')
        self.__HDUList = fits.open(filePath)
        self.__HDUList.info()
        logger.info('PRIMARY HDU header')
        print(repr(self.primary().header))
        logger.info('SPECRESP HDU header')
        print(repr(self.specresp().header))
        logger.info('SPECRESP HDU data')
        print(repr(self.specresp().data))
        print(self.specresp().data.shape)
        print(self.specresp().columns.info())

    def primary(self):
        """
        """
        return self.__HDUList['PRIMARY']

    def specresp(self):
        """
        """
        return self.__HDUList['SPECRESP']
        
        

 
def test():
    """ Test code.
    """
    import os
    from ximpol.__package__ import XIMPOL_IRF
    filePath = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.arf')
    f = xInputArfFile(filePath)
    


if __name__ == '__main__':
    test()
