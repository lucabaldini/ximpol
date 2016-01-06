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


"""Collection of system-related utilities.
"""


import shutil
import sys
import subprocess

from ximpol.utils.logging_ import logger


def cmd(cmd, verbose=False, log_file_path=None, log_file_mode='w',
        dry_run=False):
    """ Exec a system command.

    This uses subprocess internally and returns the subprocess status code
    (if the dry_run option is true the function will just print the command out
    through the logger and returns happily).

    By default the stdout and the stderr are redirected into subprocess pipes
    so that the output can be effectively used by the logger.
    It the log_file_path parameter is different from None the stdout is
    redirected to file instead.
    The rules are:

    * if verbose is True the command output is logged onto the terminal one \
    line at a time;
    * if the status code is zero we just aknowledge that before returning it;
    * upon error we try and log out both the error code and the error message.
    """
    logger.info('About to execute "%s"...' % cmd)
    if dry_run:
        logger.info('Just kidding (dry run).')
        return 0
    err = subprocess.PIPE
    if log_file_path is not None:
        out = open(log_file_path, log_file_mode)
    else:
        out = subprocess.PIPE
    process = subprocess.Popen(cmd, stdout = out, stderr = err, shell = True)
    errorCode = process.wait()
    if verbose:
        if log_file_path is None:
            output = process.stdout.read().strip('\n')
        else:
            output = open(log_file_path).read().strip('\n')
        for line in output.split('\n'):
            logger.info(line)
    if not errorCode:
        logger.info('Command executed with status code %d.' % errorCode)
    else:
        logger.error('Command returned status code %d.' % errorCode)
        logger.error('Full error message following...\n%s' %\
                         process.stderr.read().strip('\n'))
    return errorCode

def cleanup(dir_path):
    """ Remove all the files in a given folder.
    """
    file_path = os.path.join(dir_path, '*')
    cmd('rm -rf %s' % file_path)
