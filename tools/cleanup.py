#!/usr/bin/env python
#
# * Copyright (C) 2015, the ximpol team.
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


import glob
import os
import shutil
from ximpol.utils.logging_ import logger
from ximpol.utils.system_ import cmd
from ximpol import *


def cleanup(folder_path, patterns = ['*~', '*.pyc', '*.pyo']):
    """Cleanup a folder.
    """
    logger.info('Cleaning up folder %s...' % folder_path)
    fileList = []
    for pattern in patterns:
        fileList += glob.glob(os.path.join(folder_path, pattern))
    for filePath in fileList:
        logger.info('Removing %s...' % filePath)
        os.remove(filePath)

def cleanup_dist():
    """Cleanup the distribution folder.
    """
    if os.path.exists(XIMPOL_DIST):
        logger.info('Removing %s altogether...' % XIMPOL_DIST)
        shutil.rmtree(XIMPOL_DIST)
    filePath = os.path.join(XIMPOL_ROOT, 'MANIFEST')
    if os.path.exists(filePath):
        logger.info('Removing %s...' % filePath)
        os.remove(filePath)

def cleanup_doc():
    """Cleanup the doc folder.
    """
    cmd('cd %s; make clean' % XIMPOL_DOC)



if __name__ == '__main__':
    for folder_path in [
            XIMPOL_ROOT,
            XIMPOL_BIN,
            XIMPOL_CONFIG,
            XIMPOL_CORE,
            XIMPOL_DATA,
            XIMPOL_DETECTOR,
            XIMPOL_DIST,
            XIMPOL_DOC,
            XIMPOL_EVT,
            XIMPOL_EXAMPLES,
            XIMPOL_IRF,
            XIMPOL_NOTEBOOKS,
            XIMPOL_SRCMODEL,
            XIMPOL_TEST,
            XIMPOL_TOOLS,
            XIMPOL_UTILS
    ]:
        cleanup(folder_path)
    cleanup_dist()
    cleanup_doc()
