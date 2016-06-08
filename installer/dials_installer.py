from __future__ import division
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

    def rmdir(subdir):
      fullpath = os.path.join(directory, subdir)
      num_files, total_size = 0, 0
      if not os.path.exists(fullpath):
        print "Skipping", subdir
        return
      for dirpath, dirnames, filenames in os.walk(fullpath):
        for f in filenames:
          fp = os.path.join(dirpath, f)
          total_size += os.path.getsize(fp)
          num_files += 1
      print "Removing %s" % subdir
      print "  %s, %d files" % (humansize(total_size), num_files)
      shutil.rmtree(fullpath)
      print

    rmdir('base/lib/python2.7/site-packages/matplotlib-1.3.1-py2.7-linux-x86_64.egg/matplotlib/tests')
    # TODO: deduce path using something like
    # import matplotlib
    # import inspect
    # inspect.getsourcefile(matplotlib)
    # => '/scratch/wra62962/files/dials/base/lib/python2.7/site-packages/matplotlib-1.3.1-py2.7-linux-x86_64.egg/matplotlib/__init__.py'

    rmdir('base/lib/python2.7/test')
    rmdir('base/man')
    rmdir('base/share/gtk-doc')
    rmdir('modules/cbflib/doc')

if __name__ == "__main__":
  installer(sys.argv[1:]).install()
