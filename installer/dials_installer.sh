build/bin/libtbx.python modules/cctbx_project/libtbx/auto_build/create_installer.py \
	--dest            tmp/dials-installer-dev \
	--install_script  modules/dials/installer/dials_installer.py \
	--version         dev \
	--binary

cd tmp
tar -cvzf \
	dials-installer-dev.tar.gz \
	dials-installer-dev \
	--exclude dials-installer-dev/build \
	--exclude dials-installer-dev/modules \
	--exclude dials-installer-dev/base \
	--exclude dials-installer-dev/tmp

