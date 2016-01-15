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


def xpxspec(file_path, model_name):
    """Do a spectral fit in XSPEC
    """
    spec = xspec.Spectrum('test.pha')
    # These must be loaded automagically.
    spec.response = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.rmf')
    spec.response.arf = os.path.join(XIMPOL_IRF, 'fits', 'xipe_baseline.arf')
    spec.ignore('**-1.')

    model = xspec.Model(model_name)
    xspec.Fit.perform()
    xspec.Fit.show()

    xspec.Plot.device = '/xs'
    xspec.Plot.xAxis = 'keV'
    xspec.Plot('ldata')


if __name__=='__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the path to the input .pha file')
    parser.add_argument('-m', '--model', type=str, default='powerlaw',
                        help='the spectral model for the fit')

    args = parser.parse_args()
    startmsg()
    xpxspec(args.file_path, args.model)
