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


"""Unit test for the utils.profile module.
"""


import unittest

from ximpol.utils.profile import *


class testChrono(unittest.TestCase):

    """Unit test for the profile.xChrono class.
    """

    def test_basic(self):
        """Basic test of the interfaces.
        """
        c = xChrono()
        self.assertTrue(c() > 0)
        self.assertTrue(isinstance(str(c), str))


if __name__ == '__main__':
    unittest.main()
