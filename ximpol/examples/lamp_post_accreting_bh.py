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

from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedModulationCube
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.config.lamp_post_accreting_bh import pol_degree, pol_angle


CFG_FILE = os.path.join(XIMPOL_CONFIG, 'lamp_post_accreting_bh.py')
DURATION = 10000000.
E_BINNING = [1., 2., 3., 5., 7., 10.]

pipeline = xPipeline(clobber=False)
evt_file_path = pipeline.xpobssim(configfile=CFG_FILE, duration=DURATION)
pipeline.xpbin(evt_file_path, algorithm='CMAP')
mcube_file_path = pipeline.xpbin(evt_file_path, algorithm='MCUBE',
                                 ebinalg='LIST', ebinning=E_BINNING)

mcube = xBinnedModulationCube(mcube_file_path)
mcube.fit()

fig = plt.figure('Polarization degree')
mcube.plot_polarization_degree(show=False, label='Data')
pol_degree.plot(show=False, label='Model', linestyle='dashed')
plt.legend(bbox_to_anchor=(0.30, 0.95))
plt.axis([1, 10, None, None])
fig = plt.figure('Polarization angle')
mcube.plot_polarization_angle(show=False, label='Data')
pol_angle.plot(show=False, label='Model', linestyle='dashed')
plt.legend(bbox_to_anchor=(0.95, 0.95))
plt.axis([1, 10, None, None])
plt.show()
