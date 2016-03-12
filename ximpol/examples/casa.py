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
from ximpol import xpColor
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline
from ximpol.evt.binning import xBinnedMap, xBinnedModulationCube
from ximpol.utils.matplotlib_ import pyplot as plt
from matplotlib import rc
rc('text', usetex=True)


cfg_file_path = os.path.join(XIMPOL_CONFIG, 'casa.py')
evt_file_path = os.path.join(XIMPOL_DATA, 'casa.fits')
reg_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_scan.reg')
map_file_path = os.path.join(XIMPOL_DATA, 'casa_cmap.fits')
ebins_file_path = os.path.join(XIMPOL_EXAMPLES, 'casa_ebins.txt')
ebins_file_path = os.path.join(XIMPOL_EXAMPLES, 'casa_2ebins.txt')


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
    pipeline.xpobssim(configfile=cfg_file_path, duration=250000)

def select_and_bin():
    """
    """
    logger.info('Creating the mapcube for the entire source...')
    pipeline.xpbin(evt_file_path, algorithm='MCUBE', ebinalg='FILE',
                   ebinfile=ebins_file_path)
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
    fig_all = full_map.plot(show=False)

    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        #fig_all.show_circles(ra, dec, rad, lw=1)

        fig = full_map.plot(show=False, subplot=(1, 2, 1))
        plt.subplots_adjust(hspace=0.001)
        fig.show_circles(ra, dec, rad, lw=1)
        mcube_file_path = get_mcube_file_path(i)
        mcube = xBinnedModulationCube(mcube_file_path)
        mcube.plot(show=False, analyze=False, xsubplot=1)

        scale_x  = rad/numpy.cos(numpy.deg2rad(dec)) # This is to take into account the effect of the projection.
        scale_y  = rad

        for j,fit in enumerate(mcube.fit_results):
            angle = fit.phase
            angle_error = fit.phase_error
            degree = fit.polarization_degree

            #nangles=20
            #for t in range(nangles):
            #    dx = scale_x*numpy.cos(numpy.pi*2.0*t/float(nangles))#angle)
            #    dy = scale_y*numpy.sin(numpy.pi*2.0*t/float(nangles))#angle)
            #    fig.show_arrows(ra, dec, dx, dy, color='w', alpha=1, width=1,head_width=0, head_length=0)
            #    pass

            dx = scale_x*numpy.cos(angle)
            dy = scale_y*numpy.sin(angle)

            dx1 = scale_x*degree*numpy.cos(angle+angle_error)
            dy1 = scale_y*degree*numpy.sin(angle+angle_error)
            dx2 = scale_x*degree*numpy.cos(angle-angle_error)
            dy2 = scale_y*degree*numpy.sin(angle-angle_error)

            fig.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx, -dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, dx1, dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, dx2, dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)

            fig_all.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx, -dy, color=xpColor(j), alpha=1, linestyle='dashed', width=0.2,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, dx1, dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, dx2, dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)
            fig_all.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j), alpha=1, width=1,head_width=0, head_length=0)


        if save:
            fig.save(mcube_file_path.replace('.fits', '.png'))
            plt.clf()
        else:
            plt.show()
        pass
    fig_all.save(os.path.join(XIMPOL_DATA, 'casa_reg_all.png'))


if __name__ == '__main__':
    generate()
    select_and_bin()
    plot(True)
