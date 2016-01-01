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


from astropy.io import fits

from ximpol.fileio.xInputFileBase import xInputFileBase
from ximpol.__logging__ import logger


class xInputMrfFile(xInputFileBase):

    """
    """

    def __init__(self, filePath):
        """
        """
        xInputFileBase.__init__(self, filePath, 'mrf')
        self.__HDUList = fits.open(filePath)
        self.__HDUList.info()
        logger.info('PRIMARY HDU header')
        print(repr(self.primary().header))
        logger.info('MODFRESP HDU header')
        print(repr(self.modfresp().header))
        logger.info('MODFRESP HDU data')
        print(repr(self.modfresp().data))
        print(self.modfresp().data.shape)
        print(self.modfresp().columns.info())

    def primary(self):
        """
        """
        return self.__HDUList['PRIMARY']

    def modfresp(self):
        """
        """
        return self.__HDUList['MODFRESP']


def main():
    import os
    from ximpol import XIMPOL_IRF
    filePath = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.mrf')
    f = xInputMrfFile(filePath)


if __name__ == '__main__':
    main()
