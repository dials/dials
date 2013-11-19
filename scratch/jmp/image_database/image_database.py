from __future__ import division

from collections import namedtuple
ImageRecord = namedtuple(
  'ImageRecord', [
    'path',
    'mtime',
    'format_class',
    'beam',
    'detector',
    'goniometer',
    'scan'
    ])



class ImageDatabase(object):

  def __init__(self, filenames, db_name=None):

    if db_name == None:
      db_name = 'dxtbx_info.p'

    self._db_name = db_name
    self._database = self._load_files(filenames)
    self._save_database()

  def _load_files(self, filenames):
    from os.path import split, getmtime
    database = self._load_database(filenames)
    for f in filenames:
      sf = split(f)[1]
      if sf not in database or getmtime(f) != database[sf].mtime:
        database[sf] = self._load_file(f)
    return database

  def _load_file(self, filename):
    from os.path import getmtime
    return ImageRecord(
      path=filename,
      mtime=getmtime(filename),
      format_class=None,
      beam=None,
      detector=None,
      goniometer=None,
      scan=None)

  def _load_database(self, filenames):
    import cPickle as pickle
    from os.path import split, join, exists
    db_filenames = set(join(split(f)[0], self._db_name) for f in filenames)
    db = dict()
    for dbfn in db_filenames:
      if exists(dbfn):
        with open(dbfn, 'rb') as f:
          db.update(pickle.load(f))
    return db

  def _save_database(self):
    import cPickle as pickle
    with open(self._db_name, 'wb') as f:
      pickle.dump(self._database, f, protocol=pickle.HIGHEST_PROTOCOL)

if __name__ == '__main__':

  filenames = ['file1.cbf', 'file2.cbf']

  db = ImageDatabase(filenames)
  db.extract(template="*.cbf")
  db.extract(same=['beam', 'goniometer', 'detector'])
  db.extract(same=['detector'])
  sweeps = db.rotation_sweeps()
  sets = db.single_shot_sets()
