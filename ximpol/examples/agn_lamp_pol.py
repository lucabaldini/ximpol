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
from ximpol.irf.mrf import mdp99

i = 0
if i is 0:
    from ximpol.config.mcg_6_30_15 import pol_degree, pol_angle
elif i is 1:
    from ximpol.config.ark_120 import pol_degree, pol_angle
elif i is 2:
    from ximpol.config.ngc_1365 import pol_degree, pol_angle

NAME_LIST = ['mcg_6_30_15', 'ark_120', 'ngc_1365']
CFG_FILE = os.path.join(XIMPOL_CONFIG, NAME_LIST[i]+'.py')
TIME = 1000000.
E_BINNING = [2.,4.,8.]
SEED = 0

pipeline = xPipeline(clobber=True)
evt_file_path = pipeline.xpobssim(configfile=CFG_FILE, duration=TIME, seed=SEED)
pipeline.xpbin(evt_file_path, algorithm='CMAP')
mod_cube_file_path = pipeline.xpbin(evt_file_path, algorithm='MCUBE',
                                    ebinalg='LIST', ebinning=E_BINNING)
mod_cube = xBinnedModulationCube(mod_cube_file_path)
mod_cube.fit()

cnts_tot = 0
mu_tot = 0.
for i in range(0, len(mod_cube.emax)):
    cnts = mod_cube.counts[i]
    cnts_tot += cnts
    mu_tot += mod_cube.effective_mu[i]*cnts
    mdp = mdp99(mod_cube.effective_mu[i], cnts)
    logger.info('%.2f--%.2f keV: %d counts in %d s, mu %.3f, MDP %.2f%%' %\
                            (mod_cube.emin[i], mod_cube.emax[i], cnts, TIME,
                            mod_cube.effective_mu[i], 100*mdp))
mu_tot /= cnts_tot
mdp_tot = mdp99(mu_tot, cnts_tot)
logger.info('%.2f--%.2f keV: %d counts in %d s, mu %.3f, MDP %.2f%%' %\
        (mod_cube.emin[0], mod_cube.emax[len(mod_cube.emax)-1], cnts_tot,
        TIME, mu_tot, 100*mdp_tot))
logger.info('Done.')

fig = plt.figure('Polarization degree')
pol_deg_err = [r.polarization_degree_error for r in mod_cube.fit_results]
pol_deg = [r.polarization_degree for r in mod_cube.fit_results]
if pol_deg[0] < 3*pol_deg_err[0]:
    mod_cube.plot_polarization_degree(show=False, label='Data', color='gray')
else:
    mod_cube.plot_polarization_degree(show=False, label='Data')
pol_degree.plot(show=False, label='Model', linestyle='dashed', color='green')
plt.legend(bbox_to_anchor=(0.30, 0.95))
plt.axis([1, 10, 0, 0.1])

fig = plt.figure('Polarization angle')
if pol_deg[0] < 3*pol_deg_err[0]:
    mod_cube.plot_polarization_angle(show=False, label='Data', color='gray')
else:
    mod_cube.plot_polarization_angle(show=False, label='Data')
pol_angle.plot(show=False, label='Model', linestyle='dashed', color='green')
plt.legend(bbox_to_anchor=(0.95, 0.95))
plt.axis([1, 10, 0, 180])
plt.show()
