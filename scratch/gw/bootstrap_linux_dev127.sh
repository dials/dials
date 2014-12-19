export START=`pwd`
export HOST=http://cci.lbl.gov/dials/installers/dev-127
wget ${HOST}/dials-installer-dev-127-intel-linux-2.6-x86_64-scilinux.tar.gz
tar xvfz dials-installer-dev-127-intel-linux-2.6-x86_64-scilinux.tar.gz
cd dials-installer-dev-127
./install --prefix=${START}
cd ${START}
export BASE=${START}/dials-dev-127/base/bin
mkdir sources
mkdir build
cd sources
svn co https://svn.code.sf.net/p/cctbx/code/trunk cctbx_project
svn co https://svn.code.sf.net/p/dials/code/trunk dials
for MODULE in annlib annlib_adaptbx boost cbflib ccp4io ccp4io_adaptbx \
clipper dials gui_resources opt_resources scons tntbx; do
mv ../dials-dev-127/modules/${MODULE} .
done
cd ../build
${BASE}/python ../sources/cctbx_project/libtbx/configure.py cctbx rstbx dials spotfinder
. setpaths.sh
make
make
