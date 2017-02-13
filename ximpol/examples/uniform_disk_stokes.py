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


import os

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedStokesCube
from ximpol.utils.matplotlib_ import pyplot as plt

"""Script-wide simulation and analysis settings.
"""
CFG_FILE = os.path.join(XIMPOL_CONFIG, 'uniform_disk.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'uniform_disk')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
SIM_DURATION = 10000.
EMIN = 2.
EMAX = 8.
EBINS = 1
NXPIX = 64
NYPIX = 64
BINSZ = 10


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)

def run():
    PIPELINE.xpobssim(configfile=CFG_FILE, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)
    stokes_file = PIPELINE.xpbin(EVT_FILE_PATH, algorithm='SCUBE', nxpix=NXPIX,
                                 nypix=NYPIX, binsz=BINSZ, emin=EMIN,
                                 emax=EMAX, ebins=EBINS)
    stokes = xBinnedStokesCube(stokes_file)
    stokes.plot(ebin=0, slice=0, show=False)
    pol_list = stokes.polarization_degree_angle(ebin=0, smooth=1, sigma=3)
    fig_list = stokes.plot_polarization_degree_angle(pol_list, show=False)
    plt.show()

if __name__ == '__main__':
    run()
