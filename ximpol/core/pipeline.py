#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


"""We need to import the ximpol executables but we don't have an __init__ file
in there (and we don't want to add one).

Also, we need to import the
"""
from ximpol import XIMPOL_BIN
import sys
sys.path.append(XIMPOL_BIN)
sys.dont_write_bytecode = 1
from xpobssim import xpobssim, PARSER as XPOBSSIM_PARSER
from xpselect import xpselect, PARSER as XPSELECT_PARSER
from xpbin import xpbin, PARSER as XPBIN_PARSER
try:
    from xpxspec import xpxspec, PARSER as XPXSPEC_PARSER
except:
    print 'No XSPEC avalilable!'


from ximpol.utils.logging_ import logger
from ximpol.utils.os_ import rm


class xPipeline:

    """Class describing a simulation/analysis pipeline.

    Args
    ----
    clobber : bool or None
        Determines whether existing output files are overwritten. This global
        setting overrides whatever is passed as an argument to the single
        tools.
    """

    def __init__(self, clobber=None):
        """Constructor.
        """
        assert clobber in [None, True, False]
        self.clobber = clobber
        self.event_files = []
        self.binned_files = []

    def delete_event_files(self):
        """Remove all the event files generated by the pipeline.
        """
        for file_path in self.event_files:
            rm(file_path)

    def command_line(self, **kwargs):
        """Turn a dictionary into a string that is understood by argparse.
        """
        if self.clobber is not None:
            kwargs['clobber'] = self.clobber
        cmdline = ''
        for key, value in kwargs.items():
            # Need some extra care for lists...
            if isinstance(value, list):
                value = ('%s' % value).replace(' ', '')
            cmdline += '--%s %s ' % (key, value)
        cmdline.strip()
        return cmdline

    def xpobssim(self, **kwargs):
        """Generate an event file.

        All command-line switches accepted by xpobbsim can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split()
        kwargs = XPOBSSIM_PARSER.parse_args(switches).__dict__
        outfile = xpobssim(**kwargs)
        self.event_files.append(outfile)
        return outfile

    def xpselect(self, file_path, **kwargs):
        """Generate an event file.

        All command-line switches accepted by xpselect can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPSELECT_PARSER.parse_args(switches).__dict__
        outfile = xpselect(file_path, **kwargs)
        self.event_files.append(outfile)
        return outfile

    def xpbin(self, file_path, **kwargs):
        """Bin an event file.

        All command-line switches accepted by xpbin can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPBIN_PARSER.parse_args(switches).__dict__
        outfile = xpbin(file_path, **kwargs)
        self.binned_files.append(outfile)
        return outfile

    def xpxspec(self, file_path, **kwargs):
        """Analyze with XSPEC a PHA1 binned event file.

        All command-line switches accepted by xpxspec can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPXSPEC_PARSER.parse_args(switches).__dict__
        xpxspec(file_path, **kwargs)
