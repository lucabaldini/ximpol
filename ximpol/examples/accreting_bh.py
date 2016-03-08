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
from ximpol.utils.logging_ import logger
from ximpol.core.pipeline import xPipeline


CFG_FILE = os.path.join(XIMPOL_CONFIG, 'lamp_post_accreting_bh.py')
DURATION = 1000000.

pipeline = xPipeline()

# Generate the events.
evt_file_path = pipeline.xpobssim(configfile=CFG_FILE, duration=DURATION)

# Bin the events in different flavors.
pipeline.xpbin(evt_file_path, algorithm='CMAP')
pipeline.xpbin(evt_file_path, algorithm='PHA1')
pipeline.xpbin(evt_file_path, algorithm='LC')
pipeline.xpbin(evt_file_path, algorithm='MCUBE')
pipeline.delete_event_files()
