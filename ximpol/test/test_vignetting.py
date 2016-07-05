#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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


"""
"""


import sys
import os
import unittest

from ximpol.srcmodel.img import xFITSImage
from ximpol.irf import load_arf, DEFAULT_IRF_NAME
from ximpol.utils.logging_ import suppress_logging
from ximpol.utils.matplotlib_ import pyplot as plt


class TestVignetting(unittest.TestCase):

    """
    """

    def test_vignetting(self):
        """
        """
        from ximpol import XIMPOL_CONFIG
        file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_1p5_3p0_keV.fits')
        aeff = load_arf(DEFAULT_IRF_NAME)
        image = xFITSImage(file_path)
        



if __name__ == '__main__':
    unittest.main()
