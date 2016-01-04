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


import shutil
import os
import sys
import subprocess

from ximpol.utils.logging_ import logger



def getcmdoutput(cmd, ignoreExitCode = False):
    """ Execute a command in the shell and return its output.

    This turned out to be more tricky than expected. The commamnds module
    being deprecated in python 3.x, we want to use subprocess.
    subprocess.check_ouput() seemed like the cleaner way to achieve the
    goal, but it actually checks the exit code of the program and rises
    an exception if it's different from zero. Here we really want the
    command output (a la commands.getouput()), no matter what the exist code
    is, so we resort to intercepting both the stdout and the stderr and
    return the first of the second, depending on which is not empty.
    """
    p = subprocess.Popen(cmd, shell = True, stdout = subprocess.PIPE,
                         stderr = subprocess.PIPE)
    stdout, stderr = [item.strip('\n') for item in p.communicate()]
    if len(stderr):
        return stderr
    return stdout

def programAbsPath(programName):
    """ Return the absolute path of a shell program.

    TODO: ocurrently only implemented for GNU\Linux.
    """
    if os.name == 'posix':
        return getcmdoutput('which %s' % programName)
    else:
        return None

def programInstalled(programName):
    """ Return a tuple of two elements: a bool indicating whether a program
    is installed and its absolute path (None if the program does not exist).
    """
    absPath = programAbsPath(programName)
    exists = os.path.exists(absPath)
    if not exists:
        absPath = None
    return (exists, absPath)

def cp(source, dest, createTree = False):
    """ Copy a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to copy %s to %s...' % (source, dest))
    destFolder = os.path.dirname(dest)
    if not os.path.exists(destFolder) and createTree:
        createFolder(destFolder)
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
    """ Move a file.

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

def rm(filePath):
    """ Remove a file.

    Return 0 upon succesfull operation, 1 otherwise.
    """
    logger.info('About to remove file %s...' % filePath)
    if not os.path.exists(filePath):
        logger.info('File is not there, giving up...')
        return 0
    try:
        os.remove(filePath)
        logger.info('File succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove file (%s)' %  e)
        status = 1
    return status

def rmdir(folderPath):
    """ Remove an entire (empty or non empty) folder.
    """
    logger.info('About to remove folder %s...' % folderPath)
    try:
        shutil.rmtree(folderPath)
        logger.info('Folder succesfully removed.')
        status = 0
    except Exception as e:
        logger.error('Could not remove folder (%s)' %  e)
        status = 1
    return status

def cleanup(folderPath):
    """ Remove all the files in a given folder.
    """
    filePath = os.path.join(folderPath, '*')
    cmd('rm -rf %s' % filePath)

def createFolder(folderPath):
    """ Create a folder (unless it already exists).

    Return 0 upon succesfull operation, 1 otherwise.
    """
    if not os.path.exists(folderPath):
        logger.info('About to create folder %s...' % folderPath)
        try:
            os.makedirs(folderPath)
            logger.info('Folder succesfully created.')
            status = 0
        except Exception as e:
            logger.error('Could not create folder (%s)' % e)
            status = 1
        return status

def cmd(cmd, verbose = False, logFilePath = None, logFileMode = 'w',
        dryRun = False):
    """ Exec a command.

    This uses subprocess internally and returns the subprocess status code
    (if the dryRun option is true the function will just print the command out
    through the logger and returns happily).

    By default the stdout and the stderr are redirected into subprocess pipes
    so that the output can be effectively used by the logger. It the logFilePath
    parameter is different from None the stdout is redirected to file instead.
    The rules are:
    (*) if verbose is True the command output is logged onto the terminal one
    line at a time;
    (*) if the status code is zero we just aknowledge that before returning it;
    (*) upon error we try and log out both the error code and the error message.
    """
    logger.info('About to execute "%s"...' % cmd)
    if dryRun:
        logger.info('Just kidding (dry run).')
        return 0
    err = subprocess.PIPE
    if logFilePath is not None:
        out = open(logFilePath, logFileMode)
    else:
        out = subprocess.PIPE
    process = subprocess.Popen(cmd, stdout = out, stderr = err, shell = True)
    errorCode = process.wait()
    if verbose:
        if logFilePath is None:
            output = process.stdout.read().strip('\n')
        else:
            output = open(logFilePath).read().strip('\n')
        for line in output.split('\n'):
            logger.info(line)
    if not errorCode:
        logger.info('Command executed with status code %d.' % errorCode)
    else:
        logger.error('Command returned status code %d.' % errorCode)
        logger.error('Full error message following...\n%s' %\
                         process.stderr.read().strip('\n'))
    return errorCode



if __name__ == '__main__':
    print(createFolder('/test'))
    print(rm('/test'))
    print(cp('/test', '/tests'))
    print(cmd('ls', verbose = True))
    print(cmd('cacca'))
    print(cmd('cacca', logFilePath = 'test.log', verbose = True))
    print(rm('test.log'))
