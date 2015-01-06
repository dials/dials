# Setup installer.
# chem_data,dials_regression
build/bin/libtbx.python modules/cctbx_project/libtbx/auto_build/create_installer.py \
	--dest            tmp/dials-installer-dev-1000 \
	--install_script  modules/dials/installer/dials_installer.py \
	--pkg_dir         modules \
  --host-tag        mac-intel-osx-x86_64 \
	--version         dev-1000 \
	--binary
