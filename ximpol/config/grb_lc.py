#!/usr/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public Licensese as published by
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

import sys
import numpy
import random

from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.config.grb_swift_download import download_swift_grb_lc_file
from ximpol.config.grb_swift_download import get_all_swift_grb_names
from ximpol.config.grb_utils import parse_light_curve
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure

MIN_ENERGY = 0.3
MAX_ENERGY = 10.

#define a list with the name of the GRBs you want to plot
#If you want all the GRBs use the following line:
#grb_list = get_all_swift_grb_names()
    
def plot_swift_lc(grb_list,show=True):
    """Plots Swift GRB light curves.
    """
    plt.figure(figsize=(10, 8), dpi=80)
    plt.title('Swift XRT light curves')
    num_grb = 0
    for grb_name in grb_list:
        flux_outfile = download_swift_grb_lc_file(grb_name)
        if type(flux_outfile) is str:
            integral_flux_spline = parse_light_curve(flux_outfile)
            if integral_flux_spline != 0:
                if grb_name == 'GRB 130427A':
                    integral_flux_spline.plot(num_points=1000,logx=True,\
                                              logy=True,show=False,\
                                              color="red",linewidth=1.0)
                    num_grb += 1
                else:
                    c = random.uniform(0.4,0.8)
                    integral_flux_spline.plot(num_points=1000,logx=True,\
                                              logy=True,show=False,\
                                              color='%f'%c,linewidth=1.0)
                    num_grb += 1
    logger.info('%i GRBs included in the plot.'%num_grb)
    if show:
        plt.show()


def main(interactive=False):
    """Test the script plotting the light curve og GRB 130427A
    """
    #If you want all the GRBs use the following line:    
    #grb_list = get_all_swift_grb_names()
    grb_list = ['GRB 130427A']
    plot_swift_lc(grb_list,show=False)
    overlay_tag()
    save_current_figure('Swift_XRT_light_curves',
                        clear=False)
    if interactive:
        plt.show()

if __name__=='__main__':
    main(interactive=sys.flags.interactive)
