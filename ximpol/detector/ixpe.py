#!/usr/bin/env python
#
# Copyright (C) 2015, the ximpol team.
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


from ximpol.detector.xipe import *


def _full_path(file_name):
    return os.path.join(XIMPOL_DETECTOR, 'data', file_name)


"""IXPE-specific fields for the FITS headers. (These will be added to the
generic headers defined for the various extensions in the irf modules.)
"""
import ximpol.detector
ximpol.detector.xipe.INSTR_KEYWORDS = [
    ('TELESCOP', 'IXPE' , 'mission/satellite name'),
    ('INSTRUME', 'GPD'  , 'instrument/detector name'),
    ('DETNAM'  , 'ALL'  , 'specific detector name in use')
]


def make_all():
    """Create all the XIPE response functions.
    """
    # Effective area.
    aeff_file_path = _full_path('Area_IXPE_2_x3.asc')
    off_axis_data = [(10., _full_path('Area_IXPE_2_x3.asc'))]
    qeff_file_path = _full_path('eff_hedme8020_1atm_1cm_cuts80p_be50um_p_x.asc')
    make_arf(aeff_file_path, qeff_file_path, 'ixpe_baseline', off_axis_data)
    # Energy dispersion.
    eres_file_path = _full_path('eres_fwhm_hedme8020_1atm_1cm.asc')
    make_rmf(eres_file_path, 'ixpe_baseline')
    # Modulation factor.
    modf_file_path = _full_path('modfact_hedme8020_1atm_1cm_mng.asc')
    make_mrf(modf_file_path, 'ixpe_baseline')
    # Point-spread function.
    make_psf('ixpe_baseline')


if __name__ == '__main__':
    make_all()
