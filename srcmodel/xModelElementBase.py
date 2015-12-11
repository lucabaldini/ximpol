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



from ximpol.__logging__ import abort



class ModelElementKeyMissing(Exception):
    
    pass



class ModelElementKeyUnknown(Exception):

    pass



class ModelElementTypeError(Exception):

    pass



class xModelElementBase(dict):

    """ Small class encapsulating a dictionary providing additional capability
    for controlling the keys and their type.

    The basic idea is that we just grab Python dictionaries from a yaml file
    and provide here a mechanism to make sure that all the bits are in place
    to be able to instantiate proper class objects from those dictionaries.

    Since many model elements do have name, we provide a corresponding argument
    to the constructor, and a private class member and a method to retrieve it.
    The name defaults to None but can be used by subclasses.

    xModelElementBase has no required/optional keys and therefore corresponding
    xModelElementBase objects can only be instantiated as empty dictionaries
    (i.e., with no keyword arguments), after which they behave almost exactly as
    plain Python dictionaries.

    More interestingly, subclasses can redefine the REQUIRED_KEYS, OPTIONAL_KEYS
    and TYPE_DICT members to provide control on the data structure. 
    """

    REQUIRED_KEYS = []
    OPTIONAL_KEYS = []
    TYPE_DICT = {}

    def __init__(self, name = None, **kwargs):
        """ Constructor.
        
        Here is where all the control over the keys and values take place. 
        More specifically:
        1\ A ModelElementKeyMissing exception is thrown if any of the
        required keys is missing in the keyword arguments.
        2\ A ModelElementKeyUnknown exception is thrown if one of the keys
        in the keyword arguments is not listed in either the required or
        optional keys.
        3\ A ModelElementTypeError exception is thrown if the type of any of
        the values in the keyword arguments does not match that specified in
        the TYPE_DICT member.
        """
        self.__Name = name
        dict.__init__(self, **kwargs)
        for key in self.REQUIRED_KEYS:
            if not self.has_key(key):
                msg = 'missing key "%s" in %s.' %\
                      (key, self.__class__.__name__)
                raise ModelElementKeyMissing(msg)
        for key in self.keys():
            if key not in self.REQUIRED_KEYS + self.OPTIONAL_KEYS:
                msg = 'unknown key "%s" in %s.' %\
                      (key, self.__class__.__name__)
                raise ModelElementKeyUnknown(msg)
        for key, val in self.TYPE_DICT.items():
            if self.has_key(key) and not isinstance(self[key], val):
                msg = 'bad key type (%s) for "%s" in %s' %\
                      (type(self[key]), key, self.__class__.__name__)
                raise ModelElementTypeError(msg)

    def name(self):
        """ Return the element name.
        """
        return self.__Name

    def __getattr__(self, key):
        """ Overload __getattr__ method.

        Given an element instance of the xModelElementBase class, this is
        essentially to be able to write element.key as a shortcut for
        element['key'], which is handy at times.
        """
        return self[key]

    def __str__(self):
        """ String formatting.
        """
        _str = dict.__str__(self)
        if self.name() is not None:
            _str = '%s = %s' % (self.name(), _str)
        return _str

    

def test():
    """ Test code.
    """
    print(xModelElementBase())
    print(xModelElementBase('test'))
    try:
        print(xModelElementBase(key = 'value'))
    except ModelElementKeyUnknown as e:
        print(e)


if __name__ == '__main__':
    test()
