from __future__ import absolute_import, division

class cd(object):
  def __init__(self, path):
    self._newpath = path

  def __enter__(self):
    import os
    self._oldpath = os.getcwd()
    os.mkdir(self._newpath)
    os.chdir(self._newpath)

  def __exit__(self, etype, value, traceback):
    import os
    os.chdir(self._oldpath)


class cd_auto(cd):
  def __init__(self, filename):
    from os.path import splitext, basename
    from uuid import uuid4
    path = splitext(basename(filename))[0]
    path = "%s_%s" % (path, uuid4().hex)
    super(cd_auto, self).__init__(path)

