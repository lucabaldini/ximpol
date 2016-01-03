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


"""ximpol: an X-ray polarimetry simulation framework.
"""


import os

PACKAGE_NAME = 'ximpol'

"""Basic folder structure of the package.
"""
XIMPOL_ROOT = os.path.abspath(os.path.dirname(__file__))
XIMPOL_BASE = os.path.join(XIMPOL_ROOT, os.pardir)
XIMPOL_BIN = os.path.join(XIMPOL_ROOT, 'bin')
XIMPOL_CORE = os.path.join(XIMPOL_ROOT, 'core')
XIMPOL_DETECTOR = os.path.join(XIMPOL_ROOT, 'detector')
XIMPOL_DIST = os.path.join(XIMPOL_ROOT, 'dist')
XIMPOL_DOC = os.path.join(XIMPOL_BASE, 'doc')
XIMPOL_EVENT = os.path.join(XIMPOL_ROOT, 'event')
XIMPOL_FILEIO = os.path.join(XIMPOL_ROOT, 'fileio')
XIMPOL_IRF = os.path.join(XIMPOL_ROOT, 'irf')
XIMPOL_NOTEBOOKS = os.path.join(XIMPOL_BASE, 'notebooks')
XIMPOL_RELEASE = os.path.join(XIMPOL_BASE, 'release')
XIMPOL_SRCMODEL = os.path.join(XIMPOL_ROOT, 'srcmodel')
XIMPOL_TEST = os.path.join(XIMPOL_ROOT, 'test')
XIMPOL_UTILS = os.path.join(XIMPOL_ROOT, 'utils')


"""Version information.
"""
XIMPOL_VERSION_FILE_PATH = os.path.join(XIMPOL_ROOT, '__version__.py')

def version_info():
    """Read the tag and build date straight from the appropriate file.

    Use this when you don't want to import the module (i.e., at release time,
    when the file is changed), so that you don't have to bother with
    reloading stuff.
    """
    for line in open(XIMPOL_VERSION_FILE_PATH).readlines():
        exec(line.strip('\n'))
    return TAG, BUILD_DATE


"""Release notes file.
"""
XIMPOL_RELEASE_NOTES_PATH = os.path.join(XIMPOL_DOC, 'release_notes.rst')
