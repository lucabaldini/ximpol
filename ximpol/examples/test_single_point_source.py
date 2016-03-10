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

from ximpol import XIMPOL_CONFIG
from ximpol.core.pipeline import xPipeline


CFG_FILE = os.path.join(XIMPOL_CONFIG, 'test_single_point_source.py')
DURATION = 10000.

pipeline = xPipeline(clobber=False)
evt_file_path = pipeline.xpobssim(configfile=CFG_FILE, duration=DURATION)
ph1_file_path = pipeline.xpbin(evt_file_path, algorithm='PHA1')
pipeline.xpxspec(ph1_file_path)
