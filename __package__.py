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

""" Basic folder structure of the package.
"""
XIMPOL_ROOT = os.path.abspath(os.path.dirname(__file__))
XIMPOL_DIST = os.path.join(XIMPOL_ROOT, 'dist')
XIMPOL_DOC = os.path.join(XIMPOL_ROOT, 'doc')
XIMPOL_IRF = os.path.join(XIMPOL_ROOT, 'irf')
XIMPOL_UTILS = os.path.join(XIMPOL_ROOT, 'utils')


""" Version information.
"""
XIMPOL_VERSION_FILE_PATH = os.path.join(XIMPOL_ROOT, '__version__.py')
def versionInfo():
    """ Read the tag and build date straight from the appropriate file.
    
    Use this when you don't want to import the module (i.e., at release time,
    when the file is changed), so that you don't have to bother with
    reloading stuff.
    """
    for line in open(XIMPOL_VERSION_FILE_PATH).readlines():
        exec(line.strip('\n'))
    return TAG, BUILD_DATE


""" Release notes.
"""
XIMPOL_RELEASE_NOTES_PATH = os.path.join(XIMPOL_DOC, 'release.notes')



if __name__ == '__main__':
    print('XIMPOL_ROOT: %s' % XIMPOL_ROOT)

