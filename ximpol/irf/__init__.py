#!/usr/bin/env python
#
# Copyright (C) 2016, the ximpol team.
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


import os, sys

from ximpol import XIMPOL_IRF
from ximpol.irf.arf import xEffectiveArea
from ximpol.irf.psf import xPointSpreadFunction
from ximpol.irf.mrf import xModulationFactor
from ximpol.irf.rmf import xEnergyDispersion
from ximpol.utils.logging_ import logger
from ximpol.utils.os_ import check_input_file


def load_arf(irf_name, folder_path=None):
    """Facility to load the effective area for a given IRF set.
    """
    if folder_path is None:
        folder_path = os.path.join(XIMPOL_IRF,'fits')
    file_path = os.path.join(folder_path, '%s.arf' % irf_name)
    check_input_file(file_path, 'arf')
    return xEffectiveArea(file_path)

def load_psf(irf_name, folder_path=None):
    """Facility to load the point-spread function for a given IRF set.
    """
    if folder_path is None:
        folder_path = os.path.join(XIMPOL_IRF,'fits')
    file_path = os.path.join(folder_path, '%s.psf' % irf_name)
    check_input_file(file_path, 'psf')
    return xPointSpreadFunction(file_path)

def load_mrf(irf_name, folder_path=None):
    """Facility to load the modulation factor for a given IRF set.
    """
    if folder_path is None:
        folder_path = os.path.join(XIMPOL_IRF,'fits')
    file_path = os.path.join(folder_path, '%s.mrf' % irf_name)
    check_input_file(file_path, 'mrf')
    return xModulationFactor(file_path)

def load_rmf(irf_name, folder_path=None):
    """Facility to load the energy dispersion for a given IRF set.
    """
    if folder_path is None:
        folder_path = os.path.join(XIMPOL_IRF,'fits')
    file_path = os.path.join(folder_path, '%s.rmf' % irf_name)
    check_input_file(file_path, 'rmf')
    return xEnergyDispersion(file_path)

def load_irfs(irf_name, folder_path=None):
    """Facility to load all the instrument response functions corresponding
    to a given set.
    """
    aeff = load_arf(irf_name, folder_path)
    psf = load_psf(irf_name, folder_path)
    modf = load_mrf(irf_name, folder_path)
    edisp = load_rmf(irf_name, folder_path)
    return aeff, psf, modf, edisp
