#!/usr/bin/env python
#
# plugin.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from abc import ABCMeta, abstractmethod


class InterfaceMeta(ABCMeta):
  ''' The interface meta class. '''

  def __init__(self, name, bases, attrs):
    ''' Check each class has the name attribute. '''
    super(InterfaceMeta, self).__init__(name, bases, attrs)

    # Ensure interfaces and extensions have a name
    if 'name' not in self.__dict__:
      raise RuntimeError("%s has no member 'name'" % name)


class Interface(object):
  ''' The main interface class. '''

  __metaclass__ = InterfaceMeta

  # Dummy attribute
  name = ''

  @classmethod
  def extension(cls, name):
    ''' Get the requested extension class by name. '''
    choices = dict((ex.name, ex) for ex in cls.extensions())
    return choices[name]

  @classmethod
  def extensions(cls):
    ''' Iterate through the extensions '''

    # Check the given class is an interface
    if cls == Interface:
      raise RuntimeError('"Interface" has no extensions')
    elif Interface not in cls.__bases__:
      raise RuntimeError('%s is not an interface' % str(cls))

    # Get all the subclasses
    stack = list(cls.__subclasses__())
    while len(stack) > 0:
      cls = stack.pop()
      yield cls
      stack.extend(cls.__subclasses__())

  @staticmethod
  def interfaces():
    ''' Iterate through the interfaces. '''
    for iface in Interface.__subclasses__():
      yield iface

  @classmethod
  def phil_scope(cls):
    ''' Get the phil scope for the interface or extension. '''
    from libtbx.phil import parse
    if cls == Interface:
      raise RuntimeError('"Interface has no phil parameters"')
    elif not issubclass(cls, Interface):
      raise RuntimeError('%s is not an interface or extension' % str(cls))
    if Interface in cls.__bases__:
      master_scope = parse('%s .help=%s {}' % (cls.name, cls.__doc__))
      main_scope = master_scope.get_without_substitution(cls.name)
      assert(len(main_scope) == 1)
      main_scope = main_scope[0]
      if 'phil' in cls.__dict__:
        main_scope.adopt_scope(parse(cls.phil))
      if Interface in cls.__bases__:
        algorithm = parse('''
          algorithm = *%s
            .help = "The choice of algorithm"
            .type = choice
        ''' % ' '.join([ext.name for ext in cls.extensions()]))
        main_scope.adopt_scope(algorithm)
        for ext in cls.extensions():
          main_scope.adopt_scope(ext.phil_scope())
    else:
      if 'phil' in cls.__dict__:
        master_scope = parse('%s .help=%s {}' % (cls.name, cls.__doc__))
        main_scope = master_scope.get_without_substitution(cls.name)
        assert(len(main_scope) == 1)
        main_scope = main_scope[0]
        main_scope.adopt_scope(parse(cls.phil))
      else:
        master_scope = parse('')
    return master_scope
