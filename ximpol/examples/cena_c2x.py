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
import pyregion

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol.core.pipeline import xPipeline
from ximpol.utils.logging_ import logger
from ximpol.evt.binning import xBinnedModulationCube, xBinnedMap
from ximpol.evt.event import xEventFile
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.irf.mrf import mdp99
from ximpol.srcmodel.img import xFITSImage
from ximpol.irf import load_psf, DEFAULT_IRF_NAME

"""Script-wide simulation and analysis settings.
"""
INPUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'cena')
INPUT_FILE_PATH = '%s.fits' % INPUT_FILE_PATH_BASE
IMAGE_FITS_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'cena.fits')
#cena_jet+core.reg cointains the jet and core regions defined in the Chandra map
REG_SOURCE_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'cena_jet+core.reg')
#This is the file cointaining the famous coffin region
REG_JET_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'cena_jet2.reg')
#cena_knots.reg cointains the knots, core and ULX regions
REG_KNOTS_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'fits', 'cena_knots.reg')
#CONFIG_FILE_PATH = os.path.join(XIMPOL_CONFIG, 'm87.py')
ANALYSIS_FILE_PATH = '%s_analysis.txt' %INPUT_FILE_PATH_BASE

ACIS = 'i'
RA = 201.38912
DEC = -43.004776

ALGORITHM = ['MCUBE', 'CMAP', 'PHA1']
EMIN = 2.
EMAX = 8.
E_BINNING = [2.,4.,8.]
TIME = 1500000

IRF_NAME = 'ixpe_baseline' #DEFAULT_IRF_NAME
PSF = load_psf(IRF_NAME)

knots = pyregion.open(REG_KNOTS_FILE_PATH)
knots_name = ['Knot G', 'Knot F', 'Knot C', 'Knots A+B', 'Core', 'ULX']
jet_tot = numpy.zeros(len(E_BINNING))
bkg_tot = numpy.zeros(len(E_BINNING))

"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)

def run():
    out_file = PIPELINE.chandra2ximpol(INPUT_FILE_PATH, acis=ACIS,
            duration=TIME, regfile=REG_SOURCE_FILE_PATH, irfname=IRF_NAME)
            #mc=False, configfile=CONFIG_FILE_PATH)
    cmap_file = PIPELINE.xpbin(out_file, algorithm=ALGORITHM[1], nxpix=512,
                               nypix=512, binsz=1.25)
    analysis_file = open(ANALYSIS_FILE_PATH, 'w')
    
    #Now begin to analyze the data for each region
    for j, region in enumerate(knots):
        ra, dec, rad = knots[j].coord_list
        sel_file = out_file.replace('.fits', '_select_%i.fits' %j)
        PIPELINE.xpselect(out_file, ra=ra, dec=dec, rad=rad*60.,
                      outfile=sel_file)
        mod_file = PIPELINE.xpbin(sel_file, algorithm=ALGORITHM[0],
                                  ebinalg='LIST', ebinning=E_BINNING)
        mod_cube = xBinnedModulationCube(mod_file)
        event = xEventFile(sel_file)
        _title = 'Results for %s:' %knots_name[j]
        logger.info(_title)
        analysis_file.write(_title + '\n')
        energy = event.event_data['ENERGY']
        src_id = event.event_data['MC_SRC_ID']
        
        for i in range(0, len(mod_cube.emax)):
            #Select the data based on energy and the source id
            _energy_mask = (energy>mod_cube.emin[i])*(energy<mod_cube.emax[i])
            _jet_mask = _energy_mask*(src_id==1)
            _bkg_mask = _energy_mask*(src_id==0)
            _core_mask = _energy_mask*(src_id==2)
            cnts_jet = len(energy[_jet_mask])
            cnts_bkg = len(energy[_bkg_mask])
            cnts_core = len(energy[_core_mask])
            #print cnts_jet, cnts_bkg, cnts_core
            if j is 4:   #core
                src = cnts_core
                bkg = cnts_bkg + cnts_jet
            elif j is 5: #ULX
                src = cnts_bkg
                bkg = cnts_core + cnts_jet
            else:        #knots A+B, C, F or G
                src = cnts_jet
                bkg = cnts_bkg + cnts_core
                jet_tot[i] += src
                bkg_tot[i] += bkg
            
            #Calculate the mdp and show the results
            mdp = mdp99(mod_cube.effective_mu[i], src, bkg)
            _fmt = '%.2f--%.2f keV: %d source counts (%.1f%%) in %d s, MDP %.2f%%'
            _data = (mod_cube.emin[i], mod_cube.emax[i], src,
                    100*src/float(src+bkg), TIME, 100*mdp)
            _line = _fmt % _data
            logger.info(_line)
            analysis_file.write(_line + '\n')
        analysis_file.write('\n')
        
    #Finally analyze all the data within the jet region (knots A+B+C+F+G)
    _title = 'Results for the jet region:'
    logger.info(_title)
    analysis_file.write(_title + '\n')
    for i in range(0, len(jet_tot)):
        mdp = mdp99(mod_cube.effective_mu[i], jet_tot[i], bkg_tot[i])
        _fmt = '%.2f--%.2f keV: %d source counts (%.1f%%) in %d s, MDP %.2f%%'
        _data = (mod_cube.emin[i], mod_cube.emax[i], jet_tot[i],
                100*jet_tot[i]/float(jet_tot[i]+bkg_tot[i]), TIME, 100*mdp)
        _line = _fmt % _data
        logger.info(_line)
        analysis_file.write(_line + '\n')
    analysis_file.close()
    return cmap_file
    
def plot(cmap_file, draw_regions=True):
    #Plot first the ximpol map we got
    full_map = xBinnedMap(cmap_file)
    fig = full_map.plot(show=False)
    fig.recenter(RA, DEC, 2.6/60.)
    PSF.draw_psf_circle(fig, 0.1, 0.1, number=False)
    fig.show_colorscale(stretch='log', cmap = 'afmhot', vmin=2, vmax=100)
    fig.add_label(0.1, 0.9, 'Centaurus A (IXPE 1.5 Ms)', relative=True,
                    size='xx-large', color='white', horizontalalignment='left')
    
    #If draw_regions flag is set to True draw the regions with labels
    if draw_regions:
        for k in range (0,len(knots)):
            ra, dec, rad = knots[k].coord_list
            fig.show_circles(ra, dec, rad, color='green', lw=2)
        """
        region = pyregion.open(REG_JET_FILE_PATH)
        list_coord = region[0].coord_list
        x = numpy.array(list_coord[::2])
        y = numpy.array(list_coord[1::2])
        fig.show_polygons([numpy.array([x,y])], lw=1, color='white')
        """
        #Here the labels coordinates are hard-coded, they should be referred in
        #some way to the knots coordinates taken from the reg file
        fig.add_label(0.22,0.81,'G', relative=True, size='x-large',
                      color='white')
        fig.add_label(0.39,0.67,'F', relative=True, size='x-large',
                      color='white')
        fig.add_label(0.52,0.57,'C', relative=True, size='x-large',
                      color='white')
        fig.add_label(0.66,0.47,'A+B', relative=True, size='x-large',
                      color='white')
        fig.add_label(0.71,0.25,'Core', relative=True, size='x-large',
                      color='blue')
        fig.add_label(0.32,0.29,'ULX', relative=True, size='x-large',
                      color='white')
    
    #Plot the original Chandra map
    image = xFITSImage(IMAGE_FITS_PATH, build_cdf=False)
    fig2 = image.plot(show=False)
    fig2.recenter(RA, DEC, 2.6/60.)
    fig2.show_colorscale(stretch='log', cmap = 'afmhot', vmin=0.2, vmax=200)
    fig2.add_label(0.1, 0.9, 'Chandra 96.8 ks', relative=True,
                size='xx-large', color='white', horizontalalignment='left')
    #fig.show_contour(IMAGE_FITS_PATH, levels=[1, 2, 5, 50], colors='green',
    #                 smooth=5)
    if draw_regions:
        region = pyregion.open(REG_SOURCE_FILE_PATH)
        list_coord = region[0].coord_list
        x = numpy.array(list_coord[::2])
        y = numpy.array(list_coord[1::2])
        fig2.show_polygons([numpy.array([x,y])], color='green', lw=2)
        ra, dec, rad = region[1].coord_list
        fig2.show_circles(ra, dec, rad, color='green', lw=2)
        fig2.add_label(0.71,0.25,'Core', relative=True, size='x-large',
                      color='white')
        fig2.add_label(0.53,0.57,'Jet', relative=True, size='x-large',
                      color='white')
    plt.show()

if __name__ == '__main__':
    file = run()
    plot(file, draw_regions=True)
