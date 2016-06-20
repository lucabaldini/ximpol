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

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA
from ximpol.core.pipeline import xPipeline
from ximpol.utils.logging_ import logger
from ximpol.evt.binning import xBinnedModulationCube, xBinnedMap
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.irf.mrf import mdp99
from ximpol.srcmodel.img import xFITSImage

"""Script-wide simulation and analysis settings.
"""
INPUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'm87')
INPUT_FILE_PATH = '%s.fits' % INPUT_FILE_PATH_BASE
IMAGE_FITS_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'm87.fits')
REG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'm87_core+jet.reg')
CONFIG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'm87.py')

ACIS = 's'
RA = 187.7055
DEC = 12.392
RAD = 0.7
ALGORITHM = ['MCUBE', 'CMAP', 'PHA1']
EMIN = 2.
EMAX = 8.
EBINS = 1
TIME = 2000000
SEED = 0
PLOT = True

RA_CORE = 187.70634
DEC_CORE = 12.39109
RA_JET = 187.70233
DEC_JET = 12.392549
RAD_PSF = 11/60.


"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=True)

def plot(cmap_file):
    full_map = xBinnedMap(cmap_file)
    fig = full_map.plot(show=False)
    #fig.show_circles(RA, DEC, RAD/60., lw=1)
    fig.show_circles(RA_CORE, DEC_CORE, RAD_PSF/60., lw=1, color='white')
    fig.show_circles(RA_JET, DEC_JET, RAD_PSF/60., lw=1, color='white')
    fig.recenter(RA+0.003, DEC, 1/60.)
    fig.show_colorscale(stretch='linear', cmap = 'afmhot', vmin = 80,
                                                                    vmax=1500)
    fig.add_label(0.1, 0.9, 'XIPE 2 Ms', relative=True, size='xx-large',
                                color='white', horizontalalignment='left')
    image = xFITSImage(IMAGE_FITS_PATH, build_cdf=False)
    fig2 = image.plot(show=False)
    fig2.recenter(RA+0.003, DEC, 1/60.)
    fig2.show_colorscale(stretch='log', cmap = 'afmhot', vmin = 2, vmax=200)
    fig2.add_label(0.1, 0.9, 'Chandra 39.5 ks', relative=True, size='xx-large',
                      color='white', horizontalalignment='left')
    fig.show_contour(IMAGE_FITS_PATH, levels=[6, 10, 20, 50, 100],
                                                      colors='green', smooth=3)
def run():
    outfile = PIPELINE.chandra2ximpol(INPUT_FILE_PATH, acis=ACIS, seed = SEED,
            duration=TIME, regfile=REG_FILE_PATH, configfile=CONFIG_FILE_PATH)
    cmap_file = PIPELINE.xpbin(outfile, algorithm=ALGORITHM[1], nxpix=512,
                                                        nypix=512, binsz=1.25)
    #outfile = PIPELINE.xpselect(outfile, ra=RA, dec=DEC, rad=RAD)
    #outfile = PIPELINE.xpselect(outfile, ra=RA_CORE, dec=DEC_CORE, rad=RAD_PSF)
    outfile = PIPELINE.xpselect(outfile, ra=RA_JET, dec=DEC_JET, rad=RAD_PSF)
    mod_file = PIPELINE.xpbin(outfile, algorithm=ALGORITHM[0], emin=EMIN,
                                                        emax=EMAX, ebins=EBINS)

    mod_cube = xBinnedModulationCube(mod_file)
    mod_cube.fit()

    event = xEventFile(outfile)
    energy = event.event_data['ENERGY']
    src_id = event.event_data['MC_SRC_ID']
    jet_tot = 0
    core_tot = 0
    bkg_tot = 0
    mu_tot = 0.
    for i in range(0, len(mod_cube.emax)):
        mask_jet =\
                (energy>mod_cube.emin[i])*(energy<mod_cube.emax[i])*(src_id==1)
        mask_bkg =\
                (energy>mod_cube.emin[i])*(energy<mod_cube.emax[i])*(src_id==0)
        mask_core =\
                (energy>mod_cube.emin[i])*(energy<mod_cube.emax[i])*(src_id==2)
        cnts_jet = len(energy[mask_jet])
        cnts_bkg = len(energy[mask_bkg])
        cnts_core = len(energy[mask_core])
        jet_tot += cnts_jet
        bkg_tot += cnts_bkg
        core_tot += cnts_core
        cnts_tot = cnts_jet+cnts_bkg+cnts_core
        mu_tot += mod_cube.effective_mu[i]*cnts_tot
        mdp = mdp99(mod_cube.effective_mu[i], cnts_jet, cnts_bkg+cnts_core)
        logger.info('%.2f--%.2f keV: %d jet counts (%.1f%%) in %d s, mu %.3f, MDP %.2f%%' % (mod_cube.emin[i], mod_cube.emax[i], cnts_jet,\
                100*cnts_jet/float(cnts_tot), TIME, mod_cube.effective_mu[i],
                100*mdp))
        #print cnts_jet, cnts_bkg, cnts_core
    if EBINS > 1:
        cnts_tot = jet_tot+bkg_tot+core_tot
        mu_tot /= cnts_tot
        mdp_tot = mdp99(mu_tot, jet_tot, bkg_tot+core_tot)
        logger.info('%.2f--%.2f keV: %d jet counts (%.1f%%) in %d s, mu %.3f, MDP %.2f%%' % (mod_cube.emin[0], mod_cube.emax[len(mod_cube.emax)-1],\
            jet_tot, 100*jet_tot/float(cnts_tot), TIME, mu_tot, 100*mdp_tot))
        #print jet_tot, bkg_tot, core_tot
    if PLOT is True:
        plot(cmap_file)
        plt.show()
    logger.info('Done.')

if __name__ == '__main__':
    run()
