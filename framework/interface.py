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
    
    if 'name' not in self.__dict__:
      raise RuntimeError("%s has no member 'name'" % name)


class Interface(object):
  ''' The main interface class. '''

  __metaclass__ = InterfaceMeta

  name = ''

  @classmethod
  def create(cls, name, *args, **kwargs):
    ''' Create an extension for the interface with the given name. '''
    choices = dict((ex.name, ex) for ex in cls.extensions())
    return choices[name](*args, **kwargs)

  @classmethod
  def extensions(cls):
    ''' Iterate through the extensions '''
    if cls == Interface:
      raise RuntimeError('"Interface" has no extensions')
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
