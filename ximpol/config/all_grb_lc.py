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

import numpy
import random

from ximpol import XIMPOL_CONFIG
from ximpol.utils.logging_ import logger
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.config.grb_swift_download import download_swift_grb_lc_file
from ximpol.config.grb_swift_download import get_all_swift_grb_names
from ximpol.config.grb_swift_download import parse_light_curve

GRB_RA = 173.1362
GRB_DEC = 27.7129
MIN_ENERGY = 0.3
MAX_ENERGY = 10.
PL_INDEX = 1.66

def main():
    """
    """
    plt.figure(figsize=(10, 6), dpi=80)
    plt.title('Swift XRT light curves of GRBs up to 130427A')
    #get_all_swift_grb_names = ['GRB 130427A','GRB 041223']
    num_grb = 0
    for i,grb_name in enumerate(get_all_swift_grb_names()):
        flux_outfile = download_swift_grb_lc_file(grb_name)
        if type(flux_outfile) is str:
            integral_flux_spline = parse_light_curve(flux_outfile)
            if integral_flux_spline != 0:
                if grb_name == 'GRB 130427A':
                    integral_flux_spline.plot(num_points=1000,logx=True,\
                                              logy=True,show=False,\
                                              color="red",linewidth=1.0)
                    num_grb += 1
                    break
                    plt.title('Swift XRT light curves of GRBs up to now')
                else:
                    c = random.uniform(0.4,0.8)
                    integral_flux_spline.plot(num_points=1000,logx=True,\
                                              logy=True,show=False,\
                                              color='%f'%c,linewidth=1.0)
                    num_grb += 1
    print num_grb
    plt.show()
if __name__=='__main__':
    main()
