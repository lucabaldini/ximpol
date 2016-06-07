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


CFG_FILE = os.path.join(XIMPOL_CONFIG, 'casa.py')
DURATION = 250000.
E_BINNING = [1., 4., 8.]


evt_file_path = os.path.join(XIMPOL_DATA, 'casa.fits')
map_file_path = os.path.join(XIMPOL_DATA, 'casa_cmap.fits')


pipeline = xPipeline(clobber=False)

ra_min = 351.000
ra_max =  350.733

dec_min = 58.747
dec_max = 58.877

num_points = 5
radius = 0.5

OUTPUT_FILE = 'casa_polmap_info_Ebins_%s_%sarcmin.txt'%(num_points,radius)

_ra = numpy.linspace(ra_min, ra_max, num_points)
_dec = numpy.linspace(dec_min, dec_max, num_points)



def get_sel_file_path(i,j):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_polmap_Ebins%02d%02d.fits' % (i,j))

def get_mcube_file_path(i,j):
    """
    """
    return os.path.join(XIMPOL_DATA, 'casa_polmap_Ebins%02d%02d_mcube.fits' % (i,j))


def generate():
    """
    """
    pipeline.xpobssim(configfile=CFG_FILE, duration=DURATION,
                      outfile=evt_file_path)

def get_output_file():
    """
    """
    
    return os.path.join(XIMPOL_DATA, OUTPUT_FILE)
    
#loop over ra and dec (two loops) and run xpselect and xpbin to get yourself some small regions. The radius of the region you are selecting should be of the order of the psf so roughly 20 arcsecs. Then make the mcube for each of these regions and fetch the pol degree. At this point you will have arrays for ra and dec of 1-d (n and m) and the array for the pol degree needs to be of dimensions nxm.

def select_and_bin():
    
    for i,ra in enumerate(_ra):
        for j,dec in enumerate(_dec):
            #radius = 0.25 #roughly 12 arcseconds in arcminutes
            radius = 0.5 #roughly 30 arcseconds in arcminutes
         
            logger.info('Analyzing region at ra = %s, dec = %s' % (ra, dec))
            sel_file_path = get_sel_file_path(i,j)
            mcube_file_path = get_mcube_file_path(i,j)
            logger.info('Going to use %s and %s for the outputfiles..'%(sel_file_path,mcube_file_path))
            pipeline.xpselect(evt_file_path, ra=ra, dec=dec, rad=radius,
                              outfile=sel_file_path)
            pipeline.xpbin(sel_file_path, algorithm='MCUBE', ebinalg='LIST',
                           ebinning=E_BINNING, outfile = mcube_file_path)

def fit_pol_map():
    
    
    output_file =  file('%s/%s'%(XIMPOL_DATA,OUTPUT_FILE),'w')
    
    _ra = numpy.linspace(ra_min, ra_max, num_points)
    _dec = numpy.linspace(dec_min, dec_max, num_points)
    for i,ra in enumerate(_ra):
        for j,dec in enumerate(_dec):
            mcube_file_path = get_mcube_file_path(i,j)
            mcube = xBinnedModulationCube(mcube_file_path)
            mcube.fit()
            for j,fit in enumerate(mcube.fit_results):

                degree = fit.polarization_degree
                degree_error = fit.polarization_degree_error
                angle = fit.phase
                angle_error = fit.phase_error
                counts = mcube.counts[j]
                emin = mcube.emin[j]
                emax = mcube.emax[j]
                
                if degree>1.0:
                    degree = 0
                else:
                    degree = degree*100
                    degree_error = degree_error*100
                if counts < 10000:
                    degree = 0
                if degree > 50:
                    mcube.plot_bin(0)
              
                line = '%s,%s,%s,%s,%s,%s,%s\n'%(emin,emax,degree, degree_error, angle, angle_error, counts)
                print "Info",line
                output_file.write(line)

    output_file.close()
    


def plot_pol_map_from_ascii():
    
    output_file = get_output_file()
    logger.info('Opening file %s for plotting...' % output_file)
    emin,emax, degree, degree_error, angle, angle_error, counts = numpy.loadtxt(output_file, delimiter=',', unpack=True)
    
        
    _ra = numpy.linspace(ra_min, ra_max, num_points)
    _dec = numpy.linspace(dec_min, dec_max, num_points)

    _ra, _dec = numpy.meshgrid(_ra, _dec)
    fig = plt.figure(figsize=(15,15),facecolor='w')
    for i in range(len(E_BINNING) - 1):
        sigma_array = degree[emin==E_BINNING[i]]/degree_error[emin==E_BINNING[i]]
        sigma_array = sigma_array.reshape((num_points,num_points))
        pol_array = degree[emin==E_BINNING[i]].reshape((num_points,num_points))
        
        ax1 = fig.add_subplot(2,1,1)
        label_pol = 'Pol degree %.2f-%.2f keV' % (E_BINNING[i],E_BINNING[i+1])
        plt.contourf(_ra,_dec,pol_array)
        plt.colorbar()
        plt.text(0.02, 0.92, label_pol, transform=plt.gca().transAxes,
                 fontsize=25)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.show()
        ax2 = fig.add_subplot(2,2,1)
        label_sigma = 'Pol sigma %.2f-%.2f keV' % (E_BINNING[i],E_BINNING[i+1])
        plt.contourf(_ra,_dec,sigma_array)
        plt.colorbar()
        plt.text(0.02, 0.92, label_sigma, transform=plt.gca().transAxes,
                fontsize=25)
        plt.xlabel('RA')
        plt.ylabel('DEC')
        plt.show()
    
if __name__ == '__main__':
    
    #generate()
    select_and_bin()
    fit_pol_map()
    plot_pol_map_from_ascii()
