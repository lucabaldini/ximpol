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


__description__ = 'Perform a spectral fit in XSPEC'


import os
import xspec

from ximpol import XIMPOL_IRF
from ximpol.utils.logging_ import logger, startmsg


"""Command-line switches.
"""
import ast
import argparse

formatter = argparse.ArgumentDefaultsHelpFormatter
PARSER = argparse.ArgumentParser(description=__description__)
PARSER.add_argument('phafile', type=str,
                    help='the path to the input .pha file')
PARSER.add_argument('--model', type=str, default='powerlaw',
                    help='the spectral model for the fit')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')


def xpxspec(file_path, **kwargs):
    """Do a spectral fit in XSPEC
    """
    spec = xspec.Spectrum(file_path)
    # These must be loaded automagically.
    spec.response = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.rmf')
    spec.response.arf = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.arf')
    spec.ignore('**-0.5')
    spec.ignore('10.-**')
    model = xspec.Model(kwargs['model'])
    xspec.Fit.perform()
    xspec.Fit.show()
    xspec.Plot.device = '/xs'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot('ldata', 'resid')


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpxspec(args.phafile, **args.__dict__)
