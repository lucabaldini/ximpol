#!/usr/bin/env python
#
# * Copyright (C) 2015, the ximpol team.
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


import time
import os
import ximpol.__utils__ as __utils__

from ximpol.utils.logging_ import logger
from ximpol import XIMPOL_VERSION_FILE_PATH, version_info,\
    XIMPOL_RELEASE_NOTES_PATH, XIMPOL_DIST, XIMPOL_ROOT


BUILD_DATE = time.strftime('%a, %d %b %Y %H:%M:%S %z')
TAG_MODES = ['major', 'minor', 'patch']


def updateVersionInfo(mode, dryRun = False):
    """ Update the __tag__.py module with the new tag and build date.
    """
    prevTag, prevBuildDate = version_info()
    logger.info('Previous tag was %s...' % prevTag)
    version, release, patch = [int(item) for item in prevTag.split('.')]
    if mode == 'major':
        version += 1
        release = 0
        patch = 0
    elif mode == 'minor':
        release += 1
        patch = 0
    elif mode == 'patch':
        patch += 1
    else:
        abort('Unknown release mode %s.' % mode)
    nextTag = '%s.%s.%s' % (version, release, patch)
    logger.info('Writing new tag (%s) to %s...' %\
                (nextTag, XIMPOL_VERSION_FILE_PATH))
    if not dryRun:
        outputFile = open(XIMPOL_VERSION_FILE_PATH, 'w')
        outputFile.writelines('TAG = \'%s\'\n' % nextTag)
        outputFile.writelines('BUILD_DATE = \'%s\'\n' % BUILD_DATE)
        outputFile.close()
    logger.info('Done.')
    return nextTag

def updateReleaseNotes(tag, dryRun = False):
    """ Write the new tag and build date on top of the release notes
    (which must be kept up to date during the release process).
    """
    title = 'Release notes\n=============\n\n'
    logger.info('Reading in %s...' % XIMPOL_RELEASE_NOTES_PATH)
    notes = open(XIMPOL_RELEASE_NOTES_PATH).read().strip('\n').strip(title)
    logger.info('Writing out %s...' % XIMPOL_RELEASE_NOTES_PATH)
    if not dryRun:
        outputFile = open(XIMPOL_RELEASE_NOTES_PATH, 'w')
        outputFile.writelines(title)
        outputFile.writelines('\n*ximpol (%s) - %s*\n\n' % (tag, BUILD_DATE))
        outputFile.writelines(notes)
        outputFile.close()
    logger.info('Done.')

def tagPackage(mode, dryRun = False):
    """ Tag the package.

    This means:
    (*) hg pull/update to make sure we're not missing remote modification;
    (*) figure out the target tag and update the release.notes;
    (*) commit the modifications, tag and push.
    """
    __utils__.cmd('git pull', verbose = True, dryRun = dryRun)
    __utils__.cmd('git status', verbose = True, dryRun = dryRun)
    tag = updateVersionInfo(mode, dryRun)
    updateReleaseNotes(tag, dryRun)
    msg = 'Prepare for tag %s.' % tag
    __utils__.cmd('git commit -a -m "%s"' % msg, verbose = True,
                  dryRun = dryRun)
    __utils__.cmd('git push', verbose = True, dryRun = dryRun)
    msg = 'tagging version %s' % tag
    __utils__.cmd('git tag -a %s -m "%s"' % (tag, msg), verbose = True,
                  dryRun = dryRun)
    __utils__.cmd('git push --tags', verbose = True, dryRun = dryRun)
    __utils__.cmd('git status', verbose = True, dryRun = dryRun)

def distsrc():
    """ Create a plain source distribution.
    """
    tag, buildDate = version_info()
    logger.info('Creating plain source distribution...')
    distDir = os.path.join(XIMPOL_DIST, 'src')
    srcLogFilePath = 'src.log'
    # Create the distribution.
    __utils__.cmd('python setup.py sdist --dist-dir=%s --prune' % distDir,
                  verbose = False, logFilePath = srcLogFilePath)
    # Cleanup.
    __utils__.rm(srcLogFilePath)
    __utils__.rm(os.path.join(XIMPOL_ROOT, 'MANIFEST'))
    logger.info('Done.')



if __name__ == '__main__':
    from optparse import OptionParser
    parser = OptionParser()
    parser.add_option('-t', dest = 'tagmode', type = str, default = None,
                      help = 'The release tag mode %s.' % TAG_MODES)
    parser.add_option('-n', action = 'store_true', dest = 'dryrun',
                      help = 'Dry run (i.e. do not actually do anything).')
    parser.add_option('-s', action = 'store_true', dest = 'src',
                      help = 'Create a source distribution.')
    (opts, args) = parser.parse_args()
    if not opts.tagmode and not (opts.src):
        parser.print_help()
        parser.error('Please specify at least one valid option.')
    tag = None
    if opts.tagmode is not None:
        if opts.tagmode not in TAG_MODES:
            parser.error('Invalid tag mode %s (allowed: %s)' %\
                             (opts.tagmode, TAG_MODES))
        tagPackage(opts.tagmode, opts.dryrun)
    if opts.src and not opts.dryrun:
        distsrc()
