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

from ximpol import XIMPOL_CONFIG, XIMPOL_DATA, XIMPOL_DOC
from ximpol.core.pipeline import xPipeline
from ximpol.utils.logging_ import logger
from ximpol.config.single_point_source import PL_NORM, PL_INDEX

from ximpol.irf.mrf import xStokesAccumulator
from ximpol.evt.event import xEventFile
from astropy.io import fits
from ximpol.evt.binning import xBinnedModulationCube
import numpy


"""Script-wide simulation and analysis settings.
"""
CFG_FILE = os.path.join(XIMPOL_CONFIG, 'single_point_source.py')
OUT_FILE_PATH_BASE = os.path.join(XIMPOL_DATA, 'single_point_source')
EVT_FILE_PATH = '%s.fits' % OUT_FILE_PATH_BASE
SIM_DURATION = 10000.
OUTPUT_FOLDER = os.path.join(XIMPOL_DOC, 'figures', 'showcase')
mcbube_file_path = os.path.join(XIMPOL_DATA,'single_point_source_mcube.fits')

E_BINNING = [1., 2., 3., 4., 5., 6., 8., 10.]

"""Main pipeline object.
"""
PIPELINE = xPipeline(clobber=False)

def get_sel_file_path(i):
    selected_file_path = os.path.join(XIMPOL_DATA,'single_point_source_selected_%s.fits'%i)
    return selected_file_path



def run():
    PIPELINE.xpobssim(configfile=CFG_FILE, duration=SIM_DURATION,
                      outfile=EVT_FILE_PATH)
    pha1_file_path = PIPELINE.xpbin(EVT_FILE_PATH, algorithm='PHA1')
    
    ##testing stokes parameters first find them through mcube
    PIPELINE.xpbin(EVT_FILE_PATH, algorithm='MCUBE', ebinalg='LIST',
                   ebinning=E_BINNING,
                   outfile = mcbube_file_path)
    
    mcube_file = xBinnedModulationCube(mcbube_file_path)
    mcube_file.fit()
    str_line = 'Output from Stokes:\n'
    for i in range(0, len(E_BINNING) - 1):
        emax = E_BINNING[i+1]
        emin = E_BINNING[i]
       
        selected_file_path = get_sel_file_path(i)
        PIPELINE.xpselect(EVT_FILE_PATH, emin=emin,emax=emax,
                          outfile=selected_file_path)

        #open the selected file and pass the pe_angle to the stokes
        #accumulator. Do this for each energy bin and print out str.
        sim_file = xEventFile(selected_file_path)
        phase = sim_file.hdu_list[1].data['PE_ANGLE']
        stokes = xStokesAccumulator()
    
        from ximpol.irf import load_mrf
        modf = load_mrf(sim_file.irf_name())
        energy = sim_file.hdu_list[1].data['ENERGY']
    
        _energy = energy[(energy>emin)*(energy<emax)]
        eff_mu = modf.weighted_average(_energy)
    
        stokes.fill(phase[(energy>emin)*(energy<emax)])

        visibility, dvisibility = stokes.visibility()
        phase, dphase = numpy.degrees(stokes.phase())
        pol_frac, dpol_frac = stokes.polarization_frac(eff_mu)
        
        str_line+='Visibility = %.3f +- %.3f, phase = %.2f +- %.2f deg, polarization degree = %.3f +- %.3f\n'%(visibility, dvisibility, phase, dphase, pol_frac, dpol_frac)
     
   
    print "********"
    print str_line

    
if __name__ == '__main__':
    run()
