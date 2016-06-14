#!/urs/bin/env python
#
# Copyright (C) 2015--2016, the ximpol team.
#
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU GengReral Public License as published by
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
import os

from ximpol.srcmodel.gabs import mapped_column_density_HI,\
    xpeInterstellarAbsorptionModel
from ximpol.utils.logging_ import logger, abort
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import overlay_tag, save_current_figure



def main(interactive=False):
    """Make some ISM absorption related plots.
    """
    model = xpeInterstellarAbsorptionModel()

    plt.figure()
    model.xsection.plot(logx=True, logy=True, show=False)
    save_current_figure('gabs_xsection.png', clear=False)

    plt.figure()
    model.xsection_ecube().plot(logx=True, show=False)
    save_current_figure('gabs_xsection_ecube.png', clear=False)

    plt.figure()
    ra, dec = 10.684, 41.269
    column_density = mapped_column_density_HI(ra, dec, 'LAB')
    trans = model.transmission_factor(column_density)
    trans.plot(logx=True, show=False, label='$n_H = $%.3e' % column_density)
    plt.legend(loc='upper left')
    save_current_figure('gabs_trans_lab.png', clear=False)

    plt.figure()
    for column_density in [1.e20, 1.e21, 1.e22, 1.e23]:
        trans = model.transmission_factor(column_density)
        trans.plot(logx=True, show=False, label='$n_H = $%.1e' % column_density)
    plt.legend(loc='upper left')
    save_current_figure('gabs_trans_samples.png', clear=False)
    
    if interactive:
        plt.show()

        
if __name__ == '__main__':
    main(interactive=sys.flags.interactive)
