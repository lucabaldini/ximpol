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


""" Purpose: test the function multiplication.
"""

# Load the effective area of the XIPE optics and the quantum efficiency
# of the GPD.
#
aeffPath = os.path.join(XIMPOL_DETECTOR, 'data' , 'aeff_optics_xipe_m4_x3.asc')
effPath = os.path.join(XIMPOL_DETECTOR, 'data' ,
                       'eff_hedme8020_1atm_1cm_cuts80p_be50um_p_x.asc')
aeff = xFunction1dTxtFile(aeffPath, kind = 'linear', xmin = 1, xmax = 10)
eff = xFunction1dTxtFile(effPath, kind = 'linear', xmin = 1, xmax = 10)

# Multiply the two functions to get the actual effective area.
#
prod = aeff*eff

# Plot stuff.
#
fig = plt.figure(figsize = (15, 10), facecolor = 'w')
ax1 = fig.add_subplot(3, 1, 1)
ax1.plot(aeff.x, aeff.y, 'o', color = 'blue')
ax2 = fig.add_subplot(3, 1, 2)
plt.yscale('log')
ax2.plot(eff.x, eff.y, 'o', color = 'black')
ax3 = fig.add_subplot(3, 1, 3)
plt.yscale('log')
ax3.plot(prod.x, prod.y, 'o', color = 'red')
plt.show()



