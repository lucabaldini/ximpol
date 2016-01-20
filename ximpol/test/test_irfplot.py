#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


"""Unit test for IRF plotting.
"""

import unittest
import os
import sys

from ximpol.irf import load_irfs
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure


class TestIrfPlot(unittest.TestCase):

    """
    """

    @classmethod
    def setUpClass(self):
        """Setup.
        """
        self.irf_name = 'xipe_baseline'
        self.aeff, self.psf, self.modf, self.edisp = load_irfs(self.irf_name)
        self.interactive = sys.flags.interactive

    def test_irfplot(self):
        """Plot all the instrument response functions and save images.
        """
        self.aeff.plot(show=False)
        save_current_figure('%s_aeff.png' % self.irf_name,
                            show=self.interactive)
        self.psf.plot(show=False)
        save_current_figure('%s_psf.png' % self.irf_name,
                            show=self.interactive)
        self.modf.plot(show=False)
        save_current_figure('%s_modf.png' % self.irf_name,
                            show=self.interactive)
        self.edisp.plot(show=False)
        save_current_figure('%s_edisp.png' % self.irf_name,
                            show=self.interactive)


if __name__ == '__main__':
    unittest.main(exit=not sys.flags.interactive)
