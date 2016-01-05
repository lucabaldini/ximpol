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


from ximpol.irf.base import xColDefsBase, OGIP_HEADER_SPECS


"""Header specifications for the SPECRESP extension of .arf FITS files.
"""
SPECRESP_HEADER_SPECS = [
    ('EXTNAME' , 'SPECRESP', 'name of this binary table extension'),
    ('HDUCLAS1', 'RESPONSE', 'dataset relates to spectral response'),
    ('HDUCLAS2', 'SPECRESP', 'dataset contains spectral response')
] + OGIP_HEADER_SPECS


class xColDefsSPECRESP(xColDefsBase):

    """ximpol.irf.base.xColDefsBase subclass for the SPECRESP extension
    of .arf FITS files.
    """

    COLUMN_SPECS = [
        ('ENERG_LO', 'E', 'keV'),
        ('ENERG_HI', 'E', 'keV'),
        ('SPECRESP', 'E', 'cm**2')
    ]
