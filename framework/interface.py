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
from __future__ import absolute_import, division
from abc import ABCMeta, abstractmethod # implicit import


class InterfaceMeta(ABCMeta):
  ''' The interface meta class.

  This class adds some definition-time checking to the Interface base class to
  make sure that interfaces have the required fields.

  '''

  def __init__(self, name, bases, attrs):
    ''' Check each class has the name attribute. '''
    super(InterfaceMeta, self).__init__(name, bases, attrs)

    # Ensure interfaces and extensions have a name
    if 'name' not in self.__dict__:
      raise RuntimeError("%s has no member 'name'" % name)


class Interface(object):
  ''' The interface base class.

  Interfaces can be defined for automatic registration by inheriting from this
  class.

  '''

  __metaclass__ = InterfaceMeta

  # Dummy attribute
  name = ''

  @classmethod
  def extension(cls, name):
    ''' Get the requested extension class by name.

    :param name: The name of the extension
    :returns: The extension class

    '''
    choices = dict((ex.name, ex) for ex in cls.extensions())
    return choices[name]

  @classmethod
  def extensions(cls):
    ''' Iterate through the extensions

    :returns: An iterator which loops through the list of extensions.

    '''

    # Check the given class is an interface
    if cls == Interface:
      raise RuntimeError('"Interface" has no extensions')
    elif Interface not in cls.__bases__:
      raise RuntimeError('%s is not an interface' % str(cls))

    # Get all the subclasses
    stack = list(cls.__subclasses__())
    while stack:
      cls = stack.pop()
      yield cls
      stack.extend(cls.__subclasses__())

  @staticmethod
  def interfaces():
    ''' Iterate through the interfaces.

    :returns: An iterator which loops through all the defined interfaces.

    '''
    for iface in Interface.__subclasses__():
      yield iface

  @classmethod
  def phil_scope(cls):
    ''' Get the phil scope for the interface or extension.

    :returns: The phil scope for the interface or extension

    '''
    from libtbx.phil import parse
    if cls == Interface:
      raise RuntimeError('"Interface has no phil parameters"')
    elif not issubclass(cls, Interface):
      raise RuntimeError('%s is not an interface or extension' % str(cls))
    if Interface in cls.__bases__:
      doc = '\n'.join('"%s"' % d for d in cls.__doc__)
      master_scope = parse('%s .help=%s {}' % (cls.name, doc))
      main_scope = master_scope.get_without_substitution(cls.name)
      assert(len(main_scope) == 1)
      main_scope = main_scope[0]
      if 'phil' in cls.__dict__:
        main_scope.adopt_scope(cls.phil())
      if Interface in cls.__bases__:
        def ext_names(extensions):
          names = []
          default_index = -1
          for ext in extensions:
            if 'default' in ext.__dict__:
              default_index = len(names)
            names.append(ext.name)
          if default_index < 0:
            default_index = 0
          if names:
            names[default_index] = '*' + names[default_index]
          return names
        exts = list(cls.extensions())
        if exts:
          algorithm = parse('''
            algorithm = %s
              .help = "The choice of algorithm"
              .type = choice
          ''' % ' '.join(ext_names(exts)))
          main_scope.adopt_scope(algorithm)
          for ext in exts:
            main_scope.adopt_scope(ext.phil_scope())
    else:
      if 'phil' in cls.__dict__:
        help_str = '\n'.join(['"%s"' % line for line in cls.__doc__.split()])
        master_scope = parse('%s .help=%s {}' % (cls.name, help_str))
        main_scope = master_scope.get_without_substitution(cls.name)
        assert(len(main_scope) == 1)
        main_scope = main_scope[0]
        main_scope.adopt_scope(cls.phil())
      else:
        master_scope = parse('')
    return master_scope
