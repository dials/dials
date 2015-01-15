#!/bin/sh
# Debuging tool for building installers.
build/bin/libtbx.python modules/cctbx_project/libtbx/auto_build/create_installer.py \
  --install_script  modules/dials/installer/dials_installer.py \
  --version         dev \
  --binary \
  tmp/dials-installer-dev
