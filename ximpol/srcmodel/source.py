#!/urs/bin/env python
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


import numpy

from ximpol.srcmodel.img import xFitsImage


class xSourceBase:

    """
    """

    def __init__(self, name):
        """
        """
        self.name = name

    def rvs_sky_coordinates(self, size=1):
        """
        """
        pass



class xPointSource(xSourceBase):

    """
    """

    def __init__(self, name, ra, dec):
        """
        """
        xSourceBase.__init__(self, name)
        self.ra = ra
        self.dec = dec

    def rvs_sky_coordinates(self, size=1):
        """
        """
        _ra = numpy.zeros(size)
        _ra.fill(self.ra)
        _dec = numpy.zeros(size)
        _dec.fill(self.dec)
        return (_ra, dec)


class xExtendedSource(xSourceBase):

    """
    """

    def __init__(self, name, file_path):
        """
        """
        xSourceBase.__init__(self, name)
        self.image = xFitsImage(file_path)

    def rvs_sky_coordinates(self, size=1):
        """
        """
        return self.image.rvs_coordinates(size)
