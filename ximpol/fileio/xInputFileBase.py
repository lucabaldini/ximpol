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


import os

from ximpol.__logging__ import logger, abort


class xInputFileBase:

    """Small base class for a generic input file.
    """

    def __init__(self, filePath, extension = None):
        """ Constructor.
        """
        logger.info('Opening input file %s...' % filePath)
        if extension is not None and not filePath.endswith('.%s' % extension):
            abort('This does not seem to be a .%s file' % extension)
        if not os.path.exists(filePath):
            abort('File does not exist')
        if not os.path.isfile(filePath):
            abort('File is not a file :-)')


def main():
    from ximpol import XIMPOL_DETECTOR
    filePath = os.path.join(XIMPOL_DETECTOR, 'data',
                            'aeff_optics_xipe_m4_x3.asc')
    f = xInputFileBase(filePath)
    f = xInputFileBase(filePath, 'asc')
    f = xInputFileBase(filePath, 'fits')
    f = xInputFileBase(filePath.strip('c'))
    f = xInputFileBase(XIMPOL_DETECTOR)


if __name__ == '__main__':
    main()
