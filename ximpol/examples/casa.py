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
from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.irf import load_psf, DEFAULT_IRF_NAME

CFG_FILE = os.path.join(XIMPOL_CONFIG, 'casa.py')
DURATION = 150000.
E_BINNING = [2., 4., 8.]
IRF_NAME = #DEFAULT_IRF_NAME
PSF = load_psf(IRF_NAME)

evt_file_path = os.path.join(XIMPOL_DATA, 'casa.fits')
#casa_scan.reg contains 42 regions without the global one.
#reg_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_scan.reg')
#casa_multiple.reg contains 13 regions, including the global one
#reg_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_multiple.reg')
#casa_selected.reg cointains only 2 small regions and the global one
reg_file_path = os.path.join(XIMPOL_CONFIG, 'fits', 'casa_selected.reg')
map_file_path = os.path.join(XIMPOL_DATA, 'casa_cmap.fits')


pipeline = xPipeline(clobber=False)


def get_sel_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_reg%04d.fits' % i)

def get_mcube_file_path(i):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_reg%04d_mcube.fits' % i)

def generate(vignetting=True):
    """
    """
    pipeline.xpobssim(configfile=CFG_FILE, duration=DURATION,
                      outfile=evt_file_path, irfname=IRF_NAME,
                      vignetting=vignetting)

def select_and_bin():
    """
    """
    logger.info('Creating the mapcube for the entire source...')
    pipeline.xpbin(evt_file_path, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING)
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
        pipeline.xpbin(sel_file_path, algorithm='MCUBE', ebinalg='LIST',
                       ebinning=E_BINNING, outfile = mcube_file_path)

def plot(save=False, arrows=False):
    logger.info('Plotting stuff...')
    pipeline.xpbin(evt_file_path, algorithm='CMAP', outfile=map_file_path)
    regions = pyregion.open(reg_file_path)
    full_map = xBinnedMap(map_file_path)
    fig_all = full_map.plot(show=False)
    fig_all.recenter(350.86125, 58.8175, 3.3/60.)
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        fig_all.show_circles(ra, dec, rad, lw=2, color='#7ec0ee')
    PSF.draw_psf_circle(fig_all,0.9,0.85, number=False)
    fig_all.add_label(0.05,0.95, 'XIPE %s ks'%DURATION/1000., relative=True, size='x-large',
                      color='white', horizontalalignment='left')
    if save:
        fig_all.save(os.path.join(XIMPOL_DATA, 'casa_reg_all.png'))
        plt.clf()
    else:
        plt.show()
    """
    #First make the global plot
    fig_global = full_map.plot(show=False, subplot=(1, 2, 1))
    plt.subplots_adjust(hspace=0.001)
    global_mcube_file_path = os.path.join(XIMPOL_DATA, 'casa_mcube.fits' )
    mcube = xBinnedModulationCube(global_mcube_file_path)
    mcube.plot(show=False, analyze=False, xsubplot=1, full_range=True,
               simple_stat=True)
    fig_global.save(os.path.join(XIMPOL_DATA, 'casa_global.png'))
    plt.clf()
    """
    for i, region in enumerate(regions):
        ra, dec, rad = region.coord_list
        fig = full_map.plot(show=False, subplot=(1, 2, 1))
        fig.show_circles(ra, dec, rad, lw=1.5, color='#7ec0ee')
        fig.recenter(350.86125, 58.8175, 3.3/60.)
        PSF.draw_psf_circle(fig,0.9,0.85, number=False)
        fig.add_label(0.05,0.95, 'XIPE %s ks'%DURATION/1000., relative=True, size='x-large',
                      color='white', horizontalalignment='left')
        plt.subplots_adjust(hspace=0.001)
        
        mcube_file_path = get_mcube_file_path(i)
        mcube = xBinnedModulationCube(mcube_file_path)
        #set full_range to False to hide the full energy range plot
        #set simple_stat to False to show also the advanced results
        mcube.plot(show=False, analyze=False, xsubplot=1, full_range=True,
                   simple_stat=True)
        
        #if arrows flag is set to True draw the arrows for each region
        if arrows:
            #This is to take into account the effect of the projection.
            scale_x  = rad/numpy.cos(numpy.deg2rad(dec))
            scale_y  = rad
            for j, fit in enumerate(mcube.fit_results):
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

                fig.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1,
                    linestyle='dashed', width=0.2,head_width=0, head_length=0)
                fig.show_arrows(ra, dec, -dx, -dy, color=xpColor(j), alpha=1,
                    linestyle='dashed', width=0.2,head_width=0, head_length=0)
                fig.show_arrows(ra, dec, dx1, dy1, color=xpColor(j), alpha=1,
                    width=1,head_width=0, head_length=0)
                fig.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j), alpha=1,
                    width=1,head_width=0, head_length=0)
                fig.show_arrows(ra, dec, dx2, dy2, color=xpColor(j), alpha=1,
                    width=1,head_width=0, head_length=0)
                fig.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j), alpha=1,
                    width=1,head_width=0, head_length=0)

                fig_all.show_arrows(ra, dec, dx, dy, color=xpColor(j), alpha=1,
                    linestyle='dashed', width=0.2,head_width=0, head_length=0)
                fig_all.show_arrows(ra, dec, -dx, -dy, color=xpColor(j),
                    alpha=1, linestyle='dashed', width=0.2,head_width=0,
                    head_length=0)
                fig_all.show_arrows(ra, dec, dx1, dy1, color=xpColor(j),
                    alpha=1, width=1,head_width=0, head_length=0)
                fig_all.show_arrows(ra, dec, -dx1, -dy1, color=xpColor(j),
                    alpha=1, width=1,head_width=0, head_length=0)
                fig_all.show_arrows(ra, dec, dx2, dy2, color=xpColor(j),
                    alpha=1, width=1,head_width=0, head_length=0)
                fig_all.show_arrows(ra, dec, -dx2, -dy2, color=xpColor(j),
                    alpha=1, width=1,head_width=0, head_length=0)

        if save:
            fig.save(mcube_file_path.replace('.fits', '.png'))
            plt.clf()
        else:
            plt.show()
        pass


if __name__ == '__main__':
    generate()
    select_and_bin()
    #set arrows to True to show the arrows corresponding to pol angles
    plot(save=False, arrows=True)
