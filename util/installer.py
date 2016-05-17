
# XXX This is not a standalone installer!  It must be used as part of the
# framework in libtbx/auto_build.

"""
Installer script for DIALS, based on automatically generated template.  This
must be moved to the proper location to work.
"""

from __future__ import division
import os.path
import sys
libtbx_path = os.path.join(
  os.path.abspath(os.path.dirname(os.path.dirname(__file__))), "lib")
if not libtbx_path in sys.path:
  sys.path.append(libtbx_path)
from libtbx.auto_build import install_distribution

class installer (install_distribution.installer) :
  # XXX most settings can be edited here
  product_name = "DIALS"
  dest_dir_prefix = "dials"
  make_apps = []
  configure_modules = install_distribution.installer.configure_modules + \
    ['dials', 'cbflib', 'annlib_adaptbx', 'wxtbx', "gltbx"]
  include_gui_packages = True
  remove_sources_default = False # XXX toggles removal of C++ files after building
  base_package_options = ['--dials', "--all"]
  source_packages = [ "cctbx_bundle" ] + ['dials', 'cbflib', 'annlib', 'annlib_adaptbx']
  #

  installer_dir = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

  # Various functions specific to DIALS go below here
  #
  def product_specific_finalize_install (self, log) :
    # XXX remove 'phenix' dispatchers from bin directory
    bin_dir = os.path.join(self.build_dir, "bin")
    for file_name in os.listdir(bin_dir) :
      if file_name.startswith("phenix") :
        os.remove(os.path.join(bin_dir, file_name))

if __name__ == "__main__":
  installer(sys.argv[1:])
