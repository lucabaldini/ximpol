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



import os
import numpy

from ximpol.utils.xFunction1d import xFunction1d
from ximpol.__logging__ import logger, abort



class xFunction1dTxtFile(xFunction1d):

    """ Function interpolating values read from an ascii file.
    """

    def __init__(self, filePath, kind, xmin = -numpy.inf, xmax = numpy.inf):
        """
        """
        if not os.path.exists(filePath):
            abort('Could not find input file %s' % filePath)
        if not os.path.isfile(filePath):
            abort('%s is not a file')
        logger.info('Reading data values from %s...' % filePath)
        x, y = numpy.loadtxt(filePath, unpack = True)
        xFunction1d.__init__(self, x, y, kind, xmin, xmax)



def test():
    """ Test code.
    """
    from ximpol.__package__ import XIMPOL_UTILS
    from ximpol.__utils__ import rm
    filePath = os.path.join(XIMPOL_UTILS, 'tmp.txt')
    x = numpy.linspace(0, 2*numpy.pi, 20)
    y = numpy.sin(x)
    _f = open(filePath, 'w')
    for _x, _y in zip(x, y):
        _f.write('%f\t%f\n' % (_x, _y))
    _f.close()
    f = xFunction1dTxtFile(filePath, 'linear')
    f.plot()
    rm(filePath)
    f = xFunction1dTxtFile('missing_file')


if __name__ == '__main__':
    test()
