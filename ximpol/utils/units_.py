#!/urs/bin/env python
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


FACTOR_KEV_TO_ERG = 6.2415096471204e8
FACTOR_ERGCMS_TO_MCRAB = 0.48e11


def erg2keV(val):
    """Convert erg to keV.
    """
    return val*FACTOR_KEV_TO_ERG

def keV2erg(val):
    """Convert keV to erg.
    """
    return val/FACTOR_KEV_TO_ERG

def ergcms2mcrab(val):
    """
    """
    return val*FACTOR_ERGCMS_TO_MCRAB
