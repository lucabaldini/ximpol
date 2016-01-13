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


import os

from ximpol import XIMPOL_TEST_FIGURES
from ximpol.utils.matplotlib_ import pyplot as plt


def save_current_figure(file_name, clear=True):
    """Save the current matplotlib figure in `XIMPOL_DOC_FIGURES`.

    Arguments
    ---------
    file_name : string
        The name of the output file.

    clear : bool
        If `True`, the current image is cleared after the fact.
    """
    file_path = os.path.join(XIMPOL_TEST_FIGURES, file_name)
    plt.savefig(file_path)
    if clear:
        plt.clf()
