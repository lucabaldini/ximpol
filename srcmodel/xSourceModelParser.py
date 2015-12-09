#!/usr/bin/env python
# *********************************************************************
# * Copyright (C) 2015 Luca Baldini (luca.baldini@pi.infn.it)         *
# *                                                                   *
# * For the license terms see the file LICENSE, distributed           *
# * along with this software.                                         *
# *********************************************************************
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



import yaml

from ximpol.fileio.xInputFileBase import xInputFileBase
from ximpol.__logging__ import logger
from ximpol.srcmodel.__parserlookup__ import  SOURCE_CLASS_DICT



class xSourceModelParser(xInputFileBase):

    """
    """

    def __init__(self, filePath):
        """ Constructor.
        """
        logger.info('Parsing source model configuration file %s...' % filePath)
        xInputFileBase.__init__(self, filePath, 'yaml')
        tree = yaml.load(open(filePath))
        logger.info('Raw yaml tree: %s' % tree)
        self.__SourceList = []
        for srcName, data in tree.items():
            logger.info('Parsing source "%s"...' % srcName)
            assert(len(data.items()) == 1)
            srcType, srcData = data.items()[0]
            try:
                source = SOURCE_CLASS_DICT[srcType](**srcData)
                logger.info(source)
                self.__SourceList.append(source)
            except KeyError:
                abort('Unkwnown source type "%s".' % srcType)

    def sourceList(self):
        """ Return the source list.
        """
        return self.__SourceList
    



def test():
    """ Test code.
    """
    import os
    from ximpol.__package__ import XIMPOL_SRCMODEL
    filePath = os.path.join(XIMPOL_SRCMODEL, 'yaml', 'simple_source.yaml')
    parser = xSourceModelParser(filePath)
    for i, source in enumerate(parser.sourceList()):
        logger.info('%d:, %s' % (i, source))



if __name__ == '__main__':
    test()
