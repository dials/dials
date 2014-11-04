#
# DIALS INSTALLER PACKAGING CONFIGURATION
#
# This file defines installer components specific to the DIALS distribution.
# It will be used to package the results of independently build packages (e.g.
# from BuildBot) using the common installer framework.  Runtime code is in
# cctbx_project/libtbx/auto_build/assemble_binary_installer.py.
#

from __future__ import division
from libtbx.auto_build import assemble_binary_installer
import sys

class installer_builder (assemble_binary_installer.installer_builder) :
  #---------------------------------------------------------------------
  #
  # Label for the overall package
  product_name = "DIALS"
  # prefix for installer directory and tarball
  pkg_prefix = "dials"
  # script to drive actual installer
  installer_script = "dials/util/installer.py"
  # optional license file
  license = "dials/license.txt"
  # optional directory containing additional binary packages
  bin_dir = None
  # optional text files to be placed in top-level installer directory
  readme_files = [
  ]
  # modules that we need to grab source for (in addition to the cctbx_bundle
  # package)
  source_modules = [
    "annlib",
    "annlib_adaptbx",
    "dials",
  ]
  # modules that are placed in the 'base' directory of the installer rather
  # than being bundled up with everything else.
  # XXX this is probably historical legacy rather than an actual necessity;
  # but nice to keep some flexibility in packaging
  base_modules = [
  ]
  # modules that may be part of the build, but definitely *not* desired for
  # the installer.  this can include regression directories, experimental
  # packages, or third-party code with an incompatible license.
  exclude_build_modules = [
    "dials_regression",
  ]

if (__name__ == "__main__") :
  installer_builder(sys.argv[1:])
