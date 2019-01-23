from __future__ import absolute_import, division, print_function

import os

class cd(object):
  def __init__(self, path):
    self._newpath = path

  def __enter__(self):
    self._oldpath = os.getcwd()
    os.mkdir(self._newpath)
    os.chdir(self._newpath)

  def __exit__(self, etype, value, traceback):
    os.chdir(self._oldpath)
