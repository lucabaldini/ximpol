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


import numpy
from ximpol.srcmodel.roi import xPointSource
from ximpol.srcmodel.config.crab_nebula import ROI_MODEL
from ximpol.srcmodel.spectrum import power_law, constant

pulsar = xPointSource(name='pulsar', ra=ROI_MODEL.ra, dec=ROI_MODEL.dec)
pulsar.spectrum = power_law(2., 2.)
pulsar.polarization_degree = constant(0.0)
pulsar.polarization_angle = constant(0.0)

ROI_MODEL.add_source(pulsar)
