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


"""Specifications for the FITS headers related to the OGIP standards.
(These are in common to all the fits files.)
"""
OGIP_HEADER_SPECS = [
    ('HDUCLASS', 'OGIP    ', 'format conforms to OGIP standard'),
    ('HDUVERS' , '1.1.0   ', 'Version of format (OGIP memo CAL/GEN/92-002a)'),
    ('HDUDOC'  , 'OGIP memos CAL/GEN/92-002 & 92-002a', 'Documents describing the forma'),
    ('HDUVERS1', '1.0.0   ', 'Obsolete - included for backwards compatibility'),
    ('HDUVERS2', '1.1.0   ', 'Obsolete - included for backwards compatibility')
]
