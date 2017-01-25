from __future__ import absolute_import, division
import os
import shutil
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if libtbx_path not in sys.path:
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution

class installer(install_distribution.installer):
  product_name = "DIALS"
  dest_dir_prefix = "dials"
  make_apps = []
  configure_modules = ["dials", "xia2", "iota", "prime"]
  include_gui_packages = True
  base_package_options = ["--dials"]
  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
  modules = [
    # hot
    'annlib',
    'boost',
    'scons',
    'ccp4io',
    # base
    'cbflib',
    'cctbx_project',
    'gui_resources',
    'ccp4io_adaptbx',
    'annlib_adaptbx',
    'tntbx',
    'clipper',
    # dials
    'dials',
    'xia2',
    'iota',
    'prime',
  ]
  flags = list(install_distribution.installer.flags)
  try:
    flags.remove('create_versioned_dispatchers')
  except ValueError:
    pass

  def product_specific_prepackage_hook(self, directory):
    """
    Remove irrelevant files from installer.
    """
    self.print_header('Deflating installer')

    suffixes = ['B', 'KB', 'MB', 'GB', 'TB', 'PB']
    def humansize(nbytes):
      if nbytes == 0: return '0 B'
      i = 0
      while nbytes >= 1024 and i < len(suffixes)-1:
        nbytes /= 1024.
        i += 1
      f = ('%.2f' % nbytes).rstrip('0').rstrip('.')
      return '%s %s' % (f, suffixes[i])

    self._cleaned_size, self._cleaned_files = 0, 0
    def rmdir(subdir):
      fullpath = os.path.join(directory, subdir)
      num_files, total_size = 0, 0
      if not os.path.exists(fullpath):
        print "Skipping", " " * 26, subdir
        return
      for dirpath, dirnames, filenames in os.walk(fullpath):
        for f in filenames:
          fp = os.path.join(dirpath, f)
          total_size += os.path.getsize(fp)
          num_files += 1
      print "Removing %9s, %4d files from %s" % \
          (humansize(total_size), num_files, subdir)
      shutil.rmtree(fullpath)
      self._cleaned_size = self._cleaned_size + total_size
      self._cleaned_files = self._cleaned_files + num_files

    # Deduce matplotlib path
    # (base/lib/python2.??/site-packages/matplotlib-????/matplotlib)
    # (base/Python.framework/Versions/?.?/lib/python?.?/site-packages/matplotlib-(...) on MacOS)
    try:
      import matplotlib
      import inspect
      matplotpath = os.path.dirname(os.path.dirname(inspect.getsourcefile(matplotlib)))
      relpath = []
      matplotpath, d = os.path.split(matplotpath)
      relpath.append(d)
      while d and (d != 'base'):
        matplotpath, d = os.path.split(matplotpath)
        relpath.append(d)
      if d == 'base':
        relpath.reverse()
        # delete matplotlib tests
        matplotpath = os.path.join(*relpath)
        rmdir(os.path.join(matplotpath, 'matplotlib', 'tests'))
        rmdir(os.path.join(matplotpath, 'mpl_toolkits', 'tests'))

        # ...while we're here
        sitepath = os.path.dirname(matplotpath)
        rmdir(os.path.join(sitepath, 'numpy/core/tests'))
        rmdir(os.path.join(sitepath, 'numpy/doc'))
        rmdir(os.path.join(sitepath, 'numpy/distutils/tests'))
        rmdir(os.path.join(sitepath, 'numpy/f2py/docs'))

        pythonpath = os.path.dirname(sitepath)
        rmdir(os.path.join(pythonpath, 'test'))
    except Exception:
      print "Could not deduce python package paths"

    rmdir('base/share/doc')
    rmdir('base/share/gtk-doc')
    rmdir('base/share/hdf5_examples')
    rmdir('base/share/man')
    for p in ['date_time', 'filesystem', 'program_options', 'python', 'thread']:
      rmdir(os.path.join('modules/boost/libs', p, 'example'))
      rmdir(os.path.join('modules/boost/libs', p, 'doc'))
      rmdir(os.path.join('modules/boost/libs', p, 'test'))
      rmdir(os.path.join('modules/boost/libs', p, 'tutorial'))
    rmdir('modules/cbflib/doc')
    rmdir('modules/cbflib/examples')
    rmdir('modules/cbflib/ply-3.2/doc')
    rmdir('modules/cbflib/ply-3.2/example')
    rmdir('modules/cbflib/ply-3.2/test')
    rmdir('modules/clipper/examples')
    print "-" * 60
    print "Deleted %d files, decrufting installation by %s\n" % \
        (self._cleaned_files, humansize(self._cleaned_size))

if __name__ == "__main__":
  installer(sys.argv[1:]).install()
