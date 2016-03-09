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
import pyregion
import numpy

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_EXAMPLES
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedMap, xBinnedModulationCube
from ximpol.utils.matplotlib_ import pyplot as plt


cfg_file_path = os.path.join(XIMPOL_CONFIG, 'casa.py')
evt_file_path = os.path.join(XIMPOL_DATA, 'casa.fits')
reg_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_multiple.reg')
map_file_path = os.path.join(XIMPOL_DATA, 'casa_cmap.fits')
ebins_file_path = os.path.join(XIMPOL_EXAMPLES, 'casa_ebins.txt')


pipeline = xPipeline(clobber=False)


def get_sel_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_reg%04d.fits' % i)

def get_mcube_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_reg%04d_mcube.fits' % i)

def generate():
    """
    """
    pipeline.xpobssim(configfile=cfg_file_path, duration=100000)

def select_and_bin():
    """
    """
    logger.info('Opening region file %s...' % reg_file_path)
    regions = pyregion.open(reg_file_path)
    logger.info('Found %d regions...' % len(regions))
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        rad *= 60.
        logger.info('Analyzing region at ra = %s, dec = %s' % (ra, dec))
        sel_file_path = get_sel_file_path(i)
        mcube_file_path = get_mcube_file_path(i)
        pipeline.xpselect(evt_file_path, ra=ra, dec=dec, rad=rad,
                          outfile=sel_file_path)
        pipeline.xpbin(sel_file_path, algorithm='MCUBE', ebinalg='FILE',
                       ebinfile=ebins_file_path, outfile = mcube_file_path)

def plot(save=False):
    logger.info('Plotting stuff...')
    pipeline.xpbin(evt_file_path, algorithm='CMAP', outfile=map_file_path)
    regions = pyregion.open(reg_file_path)
    full_map = xBinnedMap(map_file_path)
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        fig = full_map.plot(show=False, subplot=(1, 2, 1))
        plt.subplots_adjust(hspace=0.001)
        fig.show_circles(ra, dec, rad, lw=1)
        mcube_file_path = get_mcube_file_path(i)
        mcube = xBinnedModulationCube(mcube_file_path)
        mcube.plot(show=False, analyze=False, xsubplot=1)
        for fit in mcube.fit_results:
            angle = fit.phase
            degree = fit.polarization_degree
            scale = 10.
            dx = degree/scale*numpy.cos(angle)
            dy = degree/scale*numpy.sin(angle)
            fig.show_arrows(ra, dec, dx, dy, color='g', alpha=1, width=1)
            fig.show_arrows(ra, dec, -dx, -dy, color='g', alpha=1, width=1)
        if save:
            fig.save(mcube_file_path.replace('.fits', '.png'))
            plt.clf()
        else:
            plt.show()


if __name__ == '__main__':
    generate()
    select_and_bin()
    plot(True)
