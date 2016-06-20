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
from chandra2ximpol import chandra2ximpol, PARSER as CHANDRA2XIMPOL_PARSER
from xpbin import xpbin, PARSER as XPBIN_PARSER
from xpmdp import xpmdp, PARSER as XPMDP_PARSER

from ximpol.utils.logging_ import logger
from ximpol.utils.os_ import rm

try:
    from xpxspec import xpxspec, PARSER as XPXSPEC_PARSER
except ImportError:
    logger.warn('Could not import xpxspec.')



class xPipeline:

    """Class describing a simulation/analysis pipeline.

    Arguments
    ---------
    clobber : bool or None
        Determines whether existing output files are overwritten. This global
        setting overrides whatever is passed as an argument to the single
        tools.
    """

    def __init__(self, output_path_base=None, clobber=None):
        """Constructor.
        """
        assert clobber in [None, True, False]
        self.clobber = clobber

    def evt_file_path(self):
        """
        """
        pass

    def phase_selected_file_path(self, index, label=''):
        """
        """
        pass

    def time_selected_file_path(self, index, label=''):
        """
        """
        pass

    def pos_selected_file_path(self, index, label=''):
        """
        """
        pass

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
        return outfile

    def xpselect(self, file_path, **kwargs):
        """Generate an event file.

        All command-line switches accepted by xpselect can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPSELECT_PARSER.parse_args(switches).__dict__
        outfile = xpselect(file_path, **kwargs)
        return outfile

    def xpbin(self, file_path, **kwargs):
        """Bin an event file.

        All command-line switches accepted by xpbin can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPBIN_PARSER.parse_args(switches).__dict__
        outfile = xpbin(file_path, **kwargs)
        return outfile

    def xpxspec(self, file_path, **kwargs):
        """Analyze with XSPEC a PHA1 binned event file.

        All command-line switches accepted by xpxspec can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = XPXSPEC_PARSER.parse_args(switches).__dict__
        return xpxspec(file_path, **kwargs)
    
    def chandra2ximpol(self, file_path, **kwargs):
        """Run the Chandra-to-ximpol converter.
        
        All command-line switches accepted by chandra2ximpol can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split() + [file_path]
        kwargs = CHANDRA2XIMPOL_PARSER.parse_args(switches).__dict__
        outfile = chandra2ximpol(file_path, **kwargs)
        return outfile

    def xpmdp(self, **kwargs):
        """Calculate the MDP in energy bins for a given source model.
        
        All command-line switches accepted by xpmdp can be passed as
        keyword arguments here.
        """
        switches = self.command_line(**kwargs).split()
        kwargs = XPMDP_PARSER.parse_args(switches).__dict__
        return xpmdp(**kwargs)
