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


import numpy
import os

from ximpol.srcmodel.spectrum import int_eflux2pl_norm, power_law
from ximpol.srcmodel.spectrum import xCountSpectrum
from ximpol.irf import load_arf, load_mrf, DEFAULT_IRF_NAME
from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol import XIMPOL_SRCMODEL
from matplotlib.patches import Rectangle


E_MIN = 2.
E_MAX = 8.
PL_INDEX = 2.
FLUX_REF = 1.e-10
OBS_TIME_REF = 100.e3
OBS_PLAN_LIST_PATH = os.path.join(XIMPOL_SRCMODEL, 'ascii', 'obs_plan_list.txt')

PRIORITY_ONLY = True
GRID_COLOR = 'black'

def parse_obs_plan_list(PRIORITY_ONLY):
    """
    """
    logger.info('Parsing input file %s...' % OBS_PLAN_LIST_PATH)
    src_list = []
    input_file = open(OBS_PLAN_LIST_PATH)
    input_file.next()
    for line in input_file:
        name, line, notes = line[:15].strip(), line[15:52], line[52:]
        sed_class, flux, time, mdp = line.split()
        flux = float(flux)*1e-11
        time = float(time)
        mdp = float(mdp)
        src = {'name': name, 'flux': flux, 'time':  time, 'mdp': mdp,
               'notes': notes.replace('\n', '')}
        src_list.append(src)
    return src_list

pl_norm_ref = int_eflux2pl_norm(FLUX_REF, E_MIN, E_MAX, PL_INDEX)
logger.info('PL normalization @ %.3e erg cm^-2 s^-1: %.3e keV^-1 cm^-2 s^-1' %\
            (FLUX_REF, pl_norm_ref))
aeff = load_arf(DEFAULT_IRF_NAME)
modf = load_mrf(DEFAULT_IRF_NAME)
tsamples = numpy.array([0, OBS_TIME_REF])
ebinning = numpy.array([E_MIN, E_MAX])
energy_spectrum = power_law(pl_norm_ref, PL_INDEX)
count_spectrum = xCountSpectrum(energy_spectrum, aeff, tsamples)
mdp_table = count_spectrum.build_mdp_table(ebinning, modf)
mdp_ref = mdp_table.mdp_values()[-1]
logger.info('Reference MDP for %s s: %.3f' % (OBS_TIME_REF, mdp_ref))

obs_plan_list = parse_obs_plan_list(PRIORITY_ONLY)
mirror_list = [5, 6, 7, 8]

numpy.random.seed(1)
_color = numpy.random.random((3, len(obs_plan_list)))
numpy.random.seed(17)
_disp = numpy.random.uniform(0.9, 1.4, len(obs_plan_list))

plt.figure('Average polarization degree', (14, 10))
_x = numpy.logspace(-13, -7, 100)
for obs_time in [1.e3, 10.e3, 100.e3, 1.e6]:
    _y = 100.* mdp_ref * numpy.sqrt(OBS_TIME_REF/obs_time * FLUX_REF/_x)
    plt.plot(_x, _y, color=GRID_COLOR, ls='dashed', lw=0.5)
    _i = 53
    plt.text(_x[_i], _y[_i], '$T_{obs} =$ %d ks' % (obs_time/1000.),
             color=GRID_COLOR, rotation=-45.)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('MDP 99% CL [%]')

for j, source in enumerate(obs_plan_list):
    _x = source['flux']
    _y = source['mdp']*_disp[j]
    plt.plot(_x, _y, 'o', color=_color[:,j])
    _text = source['name']
    if source['notes'] is not '':
        _text += '\n' + source['notes']
    if j in mirror_list:
        _y *= 0.8
    else:
        _y *= 1.05
    plt.text(_x, _y, _text, color=_color[:,j],
             horizontalalignment='center', size='large')
plt.axis([1e-12, 1e-7, 0.4, 30])
#plt.savefig('obs_plan.png')

plt.show()
