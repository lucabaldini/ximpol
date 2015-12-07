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
import matplotlib.pyplot as plt

from ximpol.utils.xFunction1dTxtFile import xFunction1dTxtFile
from ximpol.__package__ import XIMPOL_DETECTOR
from ximpol.__logging__ import logger


""" Purpose: test the function optimization on the XIPE effectivea area curve.
"""


# Load the XIPE effective are in two flavors---standard and optimized.
#
filePath = os.path.join(XIMPOL_DETECTOR, 'data' , 'aeff_optics_xipe_m4_x3.asc')
stdaeff = xFunction1dTxtFile(filePath, kind = 'linear', xmin = 1)
optaeff = xFunction1dTxtFile(filePath, kind = 'linear', xmin = 1)
optaeff.optimize(rtol = 0.03, atol = 15)

# Sample the two curves over a fine grid and evaluate the maximum absolute and
# relative differences.
#
x = numpy.linspace(stdaeff.xmin(), stdaeff.xmax(), 10000)
y1 = stdaeff(x)
y2 = optaeff(x)
a = abs(y1 - y2)
amax = max(a)
axmax = x[numpy.where(a == amax)[0]]
logger.info('Maximum absolute difference: %e (@ x = %e)' % (amax, axmax))
r = 0.5*abs((y1 - y2)/(y1 + y2))
rmax = max(r)
rxmax = x[numpy.where(r == rmax)[0]]
logger.info('Maximum relative difference: %e (@ x = %e)' % (rmax, rxmax))

# Finally, plot the two.
#
plt.plot(stdaeff.x, stdaeff.y, 'o', color = 'blue')
plt.plot(optaeff.x, optaeff.y, 'o', color = 'red')
plt.show()
