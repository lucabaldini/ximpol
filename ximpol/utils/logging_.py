#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
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


"""Logging utilities, building on top of the python logging module.
"""

import logging
import sys
import time


logger = logging.getLogger('ximpol')
logger.setLevel(logging.DEBUG)


""" Configure the main terminal logger.
"""
consoleHandler = logging.StreamHandler()
consoleHandler.setLevel(logging.DEBUG)
consoleFormatter = logging.Formatter(">>> %(message)s")
consoleHandler.setFormatter(consoleFormatter)
logger.addHandler(consoleHandler)


def suppress_logging():
    """Set the logging level to critical.

    This can be imported, e.g., in the unit tests are running in background and
    don't need logging.
    """
    logger.setLevel(logging.CRITICAL)


class xFileFormatter(logging.Formatter):

    """Logging file formatter class.
    """

    def format(self, record):
        """Overloaded format method.
        """
        text = '[%s] %s' % (record.levelname, record.msg)
        return text


class xFileHandler(logging.FileHandler):

    """Logging file handler class.
    """

    def __init__(self, filePath, mode = 'a', encoding = None, delay = False):
        """Constructor.
        """
        logger.info('Opening output log file %s...' % filePath)
        logging.FileHandler.__init__(self, filePath, mode, encoding, delay)
        self.setLevel(logging.DEBUG)
        self.setFormatter(xFileFormatter())
        logger.addHandler(self)

    def close(self):
        """Close the file.
        """
        if self in logger.handlers:
            logger.removeHandler(self)
            logging.FileHandler.close(self)
            logger.info('Output log file %s closed.' % self.baseFilename)


def abort(message = ''):
    """Abort the execution (via a sys.exit) with a message.

    Use this with care, and opt for custom exceptions whenever possible.
    """
    if message != '':
        message = '%s. Abort.' % message
    else:
        message = 'Abort.'
    sys.exit(message)


def startmsg():
    """Print the start message.
    """
    from ximpol.__version__ import TAG, BUILD_DATE
    print('\n    Welcome to ximple version %s (built on %s).\n' %\
              (TAG, BUILD_DATE))
    print('    Copyright (C) 2015--2016, the ximpol team.\n\n    ximple comes with ABSOLUTELY NO WARRANTY.\n    This is free software, and you are welcome to redistribute it under certain\n    conditions. See the LICENSE file for details.\n\n    Visit https://github.com/lucabaldini/ximpol for more information.\n')



if __name__ == '__main__':
    startmsg()
