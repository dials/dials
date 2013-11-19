from __future__ import division

from collections import namedtuple
ImageRecord = namedtuple(
  'ImageRecord', [
    'filename',
    'mtime',
    'format_class',
    'beam',
    'detector',
    'goniometer',
    'scan'
    ])



class ImageDatabase(object):

  def __init__(self, filenames=None, db_name=None, check_headers=False, verbose=False):

    # True/False print out a load of info
    self._verbose = verbose

    # Set the database name
    if db_name is None:
      db_name = 'image_info.p'
    self._db_name = db_name

    # Load the database
    self._db = None
    self._load(filenames, check_headers)

  def data(self):
    from copy import deepcopy
    return deepcopy(self._db)

  def extract_sweeps(self):
    pass

  def extract_stills(self):
    pass

  def extract_all_images(self):
    pass

  def _load(self, filenames, check_headers):
    from os.path import getmtime, split
    from dxtbx.format.Registry import Registry
    from time import time

    # Get the path
    if self._verbose: print 'Extracting common path from filenames'
    if filenames is not None and len(filenames) > 0:
      self._path = self._extract_path(filenames)
    else:
      self._path = ''

    # Load the database
    if self._verbose: print 'Loading database from file'
    self._db = self._load_database(self._path)

    # Update the database
    if self._verbose: print 'Loading meta data from files'
    st = time()
    for f in filenames:
      sf = split(f)[1]

      if self._db['format_class'] == None:
        self._db['format_class'] = Registry.find(f)

      if sf not in self._db or getmtime(f) != self._db[sf].mtime:
        if self._verbose: print 'Loading file: %s' % sf
        self._db[sf] = self._load_file(self._db['format_class'], f, check_headers)
    print time() - st

    # Save the database
    if self._verbose: print 'Saving database to file'
    #self._save_database()

  def _save_database(self):
    import cPickle as pickle
    from os.path import join
    with open(self._db_name, 'wb') as f:
      pickle.dump(self._db, f, protocol=pickle.HIGHEST_PROTOCOL)

  def _load_file(self, format_class, filename, check_headers):
    from os.path import getmtime

    if check_headers:
      # Read the format and meta-data
      b, d, g, s = self._extract_metadata(format_class(filename))

    else:
      format_class = None
      b, d, g, s  = None, None, None, None

    # Create image record
    return ImageRecord(
      filename=filename,
      mtime=getmtime(filename),
      format_class=format_class,
      beam=b, detector=d, goniometer=g, scan=s)

  def _extract_metadata(self, fmt):
    try: b = fmt.get_beam()
    except Exception: b = None
    try: d = fmt.get_detector()
    except Exception: d = None
    try: g = fmt.get_goniometer()
    except Exception: g = None
    try: s = fmt.get_scan()
    except Exception: s = None
    return b, d, g, s

  def _load_database(self, path):
    import cPickle as pickle
    from os.path import split, join, exists
    db = dict()
    db['format_class'] = None
    if exists(self._db_name):
      with open(self._db_name, 'rb') as f:
        db.update(pickle.load(f))
    return db

  def _extract_path(self, filenames):
    from os.path import split, abspath
    paths = list(set(split(abspath(f))[0] for f in filenames))
    if len(paths) != 1:
      raise RuntimeError('single file directory required')
    return paths[0]


if __name__ == '__main__':

  from os.path import join
  from glob import glob
  from pprint import pprint
  import sys
  path = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'

#    filenames = ['file1.cbf', 'file2.cbf'] + glob(join(path, 'centroid_*.cbf'))
  #filenames = glob(join(path, 'centroid_*.cbf'))

  db = ImageDatabase(sys.argv[1:], check_headers=True, verbose=True)
  db.extract_sweeps()
  db.extract_stills()
  db.extract_all_images()

  #pprint(db.data(), width=80, depth=3)
