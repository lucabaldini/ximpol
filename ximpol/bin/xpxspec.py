#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
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
from ximpol.evt.fitting import xSpectralFitter
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
PARSER.add_argument('--emin', type=float, default=0.5,
                    help='minimum energy for the fit')
PARSER.add_argument('--emax', type=float, default=10.,
                    help='maximum energy for the fit')
PARSER.add_argument('--plot', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='plot the fit results')
PARSER.add_argument('--clobber', type=ast.literal_eval, choices=[True, False],
                    default=True,
                    help='overwrite or do not overwrite existing output files')


def xpxspec(file_path, **kwargs):
    """Do a spectral fit in XSPEC
    """
    fitter = xSpectralFitter(file_path, **kwargs)
    if kwargs['plot']:
        fitter.plot()
    return fitter


if __name__=='__main__':
    args = PARSER.parse_args()
    startmsg()
    xpxspec(args.phafile, **args.__dict__)
