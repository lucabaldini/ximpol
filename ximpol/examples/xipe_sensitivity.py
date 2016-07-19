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
E_MAX = 10.
PL_INDEX = 2.
FLUX_REF = 1.e-10
OBS_TIME_REF = 100.e3
BLAZAR_LIST_PATH = os.path.join(XIMPOL_SRCMODEL, 'ascii', 'blazar_list.txt')

PRIORITY_ONLY = True
GRID_COLOR = 'black'
#BLAZAR_COLOR = None

def parse_blazar_list(PRIORITY_ONLY):
    """
    """
    logger.info('Parsing input file %s...' % BLAZAR_LIST_PATH)
    src_list = []
    input_file = open(BLAZAR_LIST_PATH)
    for i in range(6):
        input_file.next()
    for line in input_file:
        if line.startswith('-'):
            if PRIORITY_ONLY is False:
                for i in range(3):
                    line = input_file.next()
                PRIORITY_ONLY = True
            else:
                return src_list
        line = line.replace('BL Lac', 'BL_Lac')
        name, line, notes = line[:12].strip(), line[12:95], line[95:]
        ra, dec, opt_class, sed_class, flux_max, flux_min, p_opt_max,\
            p_opt_min = line.split()
        flux_max = float(flux_max)
        flux_min = float(flux_min)
        p_opt_max = float(p_opt_max)
        p_opt_min = float(p_opt_min)
        if p_opt_min < 0.5:
            p_opt_min = 0.5
        if p_opt_max > 0.5:
            src = {'name': name, 'flux_min': flux_min, 'flux_max': flux_max,
                   'p_opt_max' : p_opt_max, 'p_opt_min': p_opt_min}
            src_list.append(src)


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

blazar_list = parse_blazar_list(PRIORITY_ONLY)
mirror_list = [1, 3, 9, 17, 18]

numpy.random.seed(10)
_color = numpy.random.random((3, len(blazar_list)))
#numpy.random.seed(1)
#_disp = numpy.random.uniform(0.7, 2., len(blazar_list))

plt.figure('Average polarization degree', (11, 8))
_x = numpy.logspace(-13, -7, 100)
for obs_time in [1.e3, 10.e3, 100.e3, 1.e6]:
    _y = 100.* mdp_ref * numpy.sqrt(OBS_TIME_REF/obs_time * FLUX_REF/_x)
    plt.plot(_x, _y, color=GRID_COLOR, ls='dashed', lw=0.6)
    _i = 51
    if obs_time is 1.e3:
        _x_text = _x[_i]
        _y_text = _y[_i]
    plt.text(_x_text, _y_text, '$T_{obs} =$ %d ks' % (obs_time/1000.),
             color=GRID_COLOR, rotation=-43., size=14)
    _x_text /= 10
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('MDP 99% CL [%]')

for j, blazar in enumerate(blazar_list):
    _x_max = blazar['flux_max']
    _x_min = blazar['flux_min']
    _y_max = blazar['p_opt_max']
    _y_min = blazar['p_opt_min']
    plt.plot([_x_min, _x_max, _x_max], [_y_max, _y_max, _y_min],
        color=_color[:,j], lw=1.5)
    _x_text = numpy.sqrt((_x_max)*(_x_min))
    if j in mirror_list:
        _y_text = 0.88*_y_max
    else:
        _y_text = 1.02*_y_max
    plt.text(_x_text, _y_text, blazar['name'], color=_color[:,j],
             horizontalalignment='center', size='large')
plt.axis([1e-13, 1e-8, 0.5, 50])
plt.savefig('blazar_mdp_average.png')


"""
plt.figure('Average polarization degree', (14, 10))
_x = numpy.logspace(-13, -7, 100)
for obs_time in [1.e3, 10.e3, 100.e3, 1.e6]:
    _y = 100.* mdp_ref * numpy.sqrt(OBS_TIME_REF/obs_time * FLUX_REF/_x)
    plt.plot(_x, _y, color=GRID_COLOR, ls='dashed', lw=0.5)
    _i = 51
    #plt.text(_x[_i], _y[_i], '$T_{obs} =$ %d ks' % (obs_time/1000.),
    #         color=GRID_COLOR, rotation=-42., size=12)
    if obs_time is 1.e3:
        _x_text = _x[_i]
        _y_text = _y[_i]
    plt.text(_x_text, _y_text, '$T_{obs} =$ %d ks' % (obs_time/1000.),
             color=GRID_COLOR, rotation=-43., size=13)
    _x_text /= 10
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('MDP 99% CL [%]')

for j, blazar in enumerate(blazar_list):
    _x_max = blazar['flux_max']
    _x_min = blazar['flux_min']
    _y_max = blazar['p_opt_max']
    _y_min = blazar['p_opt_min']
    _rect = Rectangle((_x_min, _y_min), _x_max-_x_min, _y_max-_y_min,
                       facecolor=_color[:,j], alpha=0.25, edgecolor=None)
    plt.gca().add_patch(_rect)
    _x_text = numpy.sqrt((_x_max)*(_x_min))
    _y_text = numpy.sqrt((_y_max)*(_y_min))*_disp[j]
    plt.text(_x_text, _y_text, blazar['name'], color=_color[:,j],
             horizontalalignment='center', size='large')
plt.axis([1e-13, 1e-8, 0.5, 50])
plt.savefig('blazar_mdp_average.png')


plt.figure('Maximum polarization degree', (14, 10))
_x = numpy.logspace(-13, -7, 100)
for obs_time in [1.e3, 10.e3, 100.e3, 1.e6]:
    _y = 100.*mdp_ref * numpy.sqrt(OBS_TIME_REF/obs_time * FLUX_REF/_x)
    plt.plot(_x, _y, color=GRID_COLOR, ls='dashed', lw=1)
    _i = 52
    plt.text(_x[_i], _y[_i], '$T_{obs} =$ %d ks' % (obs_time/1000.),
             color=GRID_COLOR, rotation=-30.)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('MDP 99% CL [%]')
for blazar in blazar_list:
    _x = [blazar['flux_min'], blazar['flux_max']]
    _y = [blazar['p_max']*numpy.random.uniform(0.9, 1.1)]*2
    _p = plt.plot(_x, _y, 'o-', ls='solid', color=BLAZAR_COLOR)
    plt.text(numpy.sqrt(_x[0]*_x[1]), _y[0]*1.04, blazar['name'],
             color=_p[0].get_color(), horizontalalignment='center')
plt.text(0.62, 0.95, 'Maximum polarization degree', fontsize=18,
         transform=plt.gca().transAxes)
plt.axis([8e-13, 1e-9, 0.5, 50])
plt.savefig('blazar_mdp_maximum.png')


plt.figure('Average polarization degree', (14, 10))
_x = numpy.logspace(-13, -7, 100)
for obs_time in [1.e3, 10.e3, 100.e3, 1.e6]:
    _y = 100.*mdp_ref * numpy.sqrt(OBS_TIME_REF/obs_time * FLUX_REF/_x)
    plt.plot(_x, _y, color=GRID_COLOR, ls='dashed', lw=1)
    _i = 52
    plt.text(_x[_i], _y[_i], '$T_{obs} =$ %d ks' % (obs_time/1000.),
             color=GRID_COLOR, rotation=-30.)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('MDP 99% CL [%]')
for blazar in blazar_list:
    _x = [blazar['flux_min'], blazar['flux_max']]
    _y = [blazar['p_ave']*numpy.random.uniform(0.9, 1.1)]*2
    _p = plt.plot(_x, _y, 'o-', ls='solid', color=BLAZAR_COLOR)
    plt.text(numpy.sqrt(_x[0]*_x[1]), _y[0]*1.04, blazar['name'],
             color=_p[0].get_color(), horizontalalignment='center')
plt.text(0.62, 0.95, 'Average polarization degree', fontsize=18,
         transform=plt.gca().transAxes)
plt.axis([8e-13, 1e-9, 0.5, 50])
plt.savefig('blazar_mdp_average.png')


plt.figure()
_x = numpy.logspace(-12, -7, 100)
for mdp in [0.01, 0.03, 0.1, 0.3]:
    _y = OBS_TIME_REF/1000. * (mdp_ref/mdp)**2. * (FLUX_REF/_x)
    plt.plot(_x, _y, color='blue', ls='dashed', lw=1)
    _i = 20
    plt.text(_x[_i], _y[_i], 'MDP = %.0f %%' % (mdp*100.), color='blue',
             rotation=-38.)
plt.xscale('log')
plt.yscale('log')
plt.xlabel('Integral energy flux %.0f-%.0f keV [erg cm$^{-2}$ s$^{-1}$]' %\
           (E_MIN, E_MAX))
plt.ylabel('Observing time [ks]')
plt.axis([None, None, 0.1, 10000])
"""

plt.show()
