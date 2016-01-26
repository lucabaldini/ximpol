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


__description__ = 'Quick FITS image viewer'


from astropy.io import fits

from ximpol.srcmodel.img import xFITSImage
from ximpol.utils.matplotlib_ import pyplot as plt
from ximpol.utils.matplotlib_ import context_no_grids
from ximpol.utils.logging_ import logger, startmsg


def xpimgview(file_path, output_file):
    """Quick FITS image viewer.
    """
    img = xFITSImage(file_path, build_cdf=False)
    img.plot(show=False)
    if output_file is not None:
        logger.info('Saving image to %s...' % output_file)
        plt.savefig(output_file)
    plt.show()



if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(description=__description__)
    parser.add_argument('file_path', type=str,
                        help='the input FITS image file')
    parser.add_argument('-o', '--output-file', type=str, default=None,
                        help='path to the output file to save the image')
    args = parser.parse_args()
    startmsg()
    xpimgview(args.file_path, args.output_file)
