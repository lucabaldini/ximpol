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


import argparse
import matplotlib.pyplot as plt
from astropy.io import fits

from ximpol.__logging__ import logger, startmsg



parser = argparse.ArgumentParser(description = '')
parser.add_argument('filePath', type = str,
                    help = 'the input .fits file containing the image')
args = parser.parse_args()

startmsg()
data = fits.getdata(args.filePath)
plt.imshow(data, cmap = 'afmhot')
plt.colorbar()
plt.show()
