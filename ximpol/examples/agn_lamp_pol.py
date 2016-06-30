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
import numpy

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
TIME = 30000000.
E_BINNING = [2.,4.,6.,8.]
SIM_NUM = 5

pol_degree_array = numpy.empty([SIM_NUM, len(E_BINNING)])
pol_degree_error_array = numpy.empty([SIM_NUM, len(E_BINNING)])
pol_angle_array = numpy.empty([SIM_NUM, len(E_BINNING)])
pol_angle_error_array = numpy.empty([SIM_NUM, len(E_BINNING)])

pipeline = xPipeline(clobber=True)
for j in range(0, SIM_NUM):
    SEED = j
    evt_file_path = pipeline.xpobssim(configfile=CFG_FILE, duration=TIME, seed=SEED)
    #pipeline.xpbin(evt_file_path, algorithm='CMAP')
    mod_cube_file_path = pipeline.xpbin(evt_file_path, algorithm='MCUBE',
                                        ebinalg='LIST', ebinning=E_BINNING)
    mod_cube = xBinnedModulationCube(mod_cube_file_path)
    mod_cube.fit()
    pol_degree_array[j,:] = numpy.array([r.polarization_degree for
                                         r in mod_cube.fit_results])
    pol_degree_error_array[j,:] = numpy.array([r.polarization_degree_error for
                                               r in mod_cube.fit_results])
    pol_angle_array[j,:] = numpy.array([numpy.degrees(r.phase) for
                                        r in mod_cube.fit_results])
    pol_angle_error_array[j,:] = numpy.array([numpy.degrees(r.phase_error) for
                                              r in mod_cube.fit_results])

for i in range(0, len(mod_cube.emax)):
    cnts = mod_cube.counts[i]
    logger.info('%.2f--%.2f keV: %d counts in %d s, mu %.3f, MDP %.2f%%' %\
            (mod_cube.emin[i], mod_cube.emax[i], mod_cube.counts[i], TIME,
             mod_cube.effective_mu[i], 100*mod_cube.mdp99[i]))
logger.info('Done.')

pol_deg = numpy.mean(pol_degree_array, axis=0)
pol_deg_err = numpy.mean(pol_degree_error_array, axis=0)
pol_ang = numpy.mean(pol_angle_array, axis=0)
pol_ang_err = numpy.mean(pol_angle_error_array, axis=0)

fig = plt.figure('Polarization degree')
plt.xlabel('Energy [keV]')
plt.ylabel('Polarization degree')
bad = pol_deg < 3*pol_deg_err
good = numpy.invert(bad)
bad[len(E_BINNING)-1] = False
good[len(E_BINNING)-1] = False
_dx = numpy.array([mod_cube.emean - mod_cube.emin, mod_cube.emax - mod_cube.emean])
if bad.sum() > 0:
    plt.errorbar(mod_cube.emean[bad], pol_deg[bad], pol_deg_err[bad],
                            _dx.T[bad].T, fmt='o', label='Data', color='gray')
if good.sum() > 0:
    plt.errorbar(mod_cube.emean[good], pol_deg[good], pol_deg_err[good],
                            _dx.T[good].T, fmt='o', label='Data', color='blue')
pol_degree.plot(show=False, label='Model', linestyle='dashed', color='green')
plt.legend(bbox_to_anchor=(0.30, 0.95))
plt.axis([1, 10, 0, 0.1])
#plt.savefig('/home/nicco/%s_0_0_10x%dMs_polarization_degree.png' %\
#                                                (NAME_LIST[0:3], TIME/100000))

fig = plt.figure('Polarization angle')
plt.xlabel('Energy [keV]')
plt.ylabel('Polarization angle [$^\circ$]')
if bad.sum() > 0:
    plt.errorbar(mod_cube.emean[bad], pol_ang[bad], pol_ang_err[bad],
                            _dx.T[bad].T, fmt='o', label='Data', color='gray')
if good.sum() > 0:
    plt.errorbar(mod_cube.emean[good], pol_ang[good], pol_ang_err[good],
                            _dx.T[good].T, fmt='o', label='Data', color='blue')
pol_angle.plot(show=False, label='Model', linestyle='dashed', color='green')
plt.legend(bbox_to_anchor=(0.95, 0.95))
plt.axis([1, 10, 0, 180])
#plt.savefig('/home/nicco/%s_0_0_10x%dMs_polarization_angle.png'%\
#                                                (NAME_LIST[0:3], TIME/100000))
plt.show()
