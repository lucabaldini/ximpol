#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
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


"""Collection of os-related utilities.
"""


import os
import shutil

from ximpol.utils.logging_ import logger


def mkdir(dir_path):
    """ Create a directory (unless it already exists).

    Return 0 upon succesfull operation, 1 otherwise.
    """
    if not os.path.exists(dir_path):
        logger.info('About to create folder %s...' % dir_path)
        try:
            os.makedirs(dir_path)
            logger.info('Folder succesfully created.')
            status = 0
        except Exception as e:
            logger.error('Could not create folder (%s)' % e)
            status = 1
        return status

def cp(source, dest, create_tree = False):
    """Copy a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to copy %s to %s...' % (source, dest))
    destFolder = os.path.dirname(dest)
    if not os.path.exists(destFolder) and createTree:
        mkdir(destFolder)
    try:
        if os.path.isdir(source):
            shutil.copytree(source, dest)
        else:
            shutil.copy(source, dest)
        logger.info('File succesfully copied.')
        status = 0
    except Exception as e:
        logger.error('Could not copy file (%s)' % e)
        status = 1
    return status

def mv(source, dest):
    """Move a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to move %s to %s...' % (source, dest))
    try:
        shutil.move(source, dest)
        logger.info('File succesfully copied.')
        status = 0
    except Exception as e:
        logger.error('Could not move file (%s)' % e)
        status = 1
    return status

def rm(file_path):
    """ Remove a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to remove file %s...' % file_path)
    if not os.path.exists(file_path):
        logger.info('File is not there, giving up...')
        return 0
    try:
        os.remove(file_path)
        logger.info('File succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove file (%s)' %  e)
        status = 1
    return status

def rmdir(dir_path):
    """ Remove an entire (empty or non empty) folder.
    """
    logger.info('About to remove folder %s...' % dir_path)
    try:
        shutil.rmtree(dir_path)
        logger.info('Folder succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove folder (%s)' %  e)
        status = 1
    return status

def cleanup(dir_path):
    """ Remove all the files in a given folder.
    """
    file_path = os.path.join(dir_path, '*')
    cmd('rm -rf %s' % file_path)
