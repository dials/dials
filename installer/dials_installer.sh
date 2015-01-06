# Setup installer.
# chem_data,dials_regression
build/bin/libtbx.python modules/cctbx_project/libtbx/auto_build/create_installer.py \
	--dest            tmp/dials-installer \
	--install_script  modules/dials/installer/dials_installer.py \
  --host-tag        mac-intel-osx-x86_64 \
	--binary
