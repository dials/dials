from __future__ import division
#!/usr/bin/env python
#
# toplevel.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from abc import ABCMeta, abstractmethod

class Interface(object):

  __metaclass__ = ABCMeta

  def __init__(self, maxiter=20):
    self._prepared = False
    self._processed = False
    self._finished = False
    self._maxiter = maxiter

  def run(self):
    count1 = 0
    while not self.finished:
      count2 = 0
      while not self.processed:
        count3 = 0
        while not self.prepared:
          self.prepare()
          count3 += 1
          if count3 >= self._maxiter:
            raise RuntimeError('maximum prepare iterations reached')
        self.process()
        count2 += 1
        if count2 >= self._maxiter:
          raise RuntimeError('maximum process iterations reached')
      self.finish()
      count1 += 1
      if count1 >= self._maxiter:
        raise RuntimeError('maximum finish iterations reached')

  @abstractmethod
  def prepare(self):
    pass

  @abstractmethod
  def process(self):
    pass

  @abstractmethod
  def finish(self):
    pass

  @property
  def prepared(self):
    return self._prepared

  @prepared.setter
  def prepared(self, value):
    self._prepared = value
    if not self.prepared:
      self.processed = False

  @property
  def processed(self):
    return self._processed

  @processed.setter
  def processed(self, value):
    self._processed = value
    if not self.processed:
      self.finished = False

  @property
  def finished(self):
    return self._finished

  @finished.setter
  def finished(self, value):
    self._finished = value
