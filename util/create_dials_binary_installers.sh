#!/bin/sh
#
# BuildBot (etc.) harvesting for DIALS binary installer.  This assumes that the
# nightly build is complete and ready for packaging.  The first argument,
# BUILD_DIR, should be a directory with the following contents:
#   base/ : Python and friends
#   build/ : compiled libraries and executables
#   modules/ : source code, etc.
#
# The intention is to make this as independent of specific setup at each site,
# but some customization is recommended, e.g. specifying the final location
# for the completed installer packages (DIALS_DEST_DIR).
#
if [ -z "${DIALS_DEST_DIR}" ]; then
  # XXX This should be defined to wherever installers live, minus version
  DIALS_DEST_DIR=""
fi
DIALS_DEBUG=1
#
# Stuff below here should be site-independent...
#
BUILD_DIR=$1
DIALS_VERSION=$2
HOST_TAG=$3
if [ -z "$BUILD_DIR" ] || [ ! -d "$BUILD_DIR" ]; then
  echo "BUILD_DIR must be first argument!"
  exit 1
fi
if [ -z "$DIALS_VERSION" ]; then
  echo "DIALS_VERSION must be second argument!"
  exit 1
fi
if [ -z "$HOST_TAG" ]; then
  echo "DIALS_VERSION must be second argument!"
  exit 1
fi
DEV_BUILD=`echo ${DIALS_VERSION} | grep -c dev`
PYTHON_BIN=${BUILD_DIR}/build/bin/libtbx.python
$PYTHON_BIN ${BUILD_DIR}/modules/dials/util/assemble_installer.py \
  --version=${DIALS_VERSION} --host-tag=${HOST_TAG} $BUILD_DIR
DIALS_INSTALLER="phenix-installer-${DIALS_VERSION}"
DIALS_TAR="phenix-installer-${DIALS_VERSION}-${HOST_TAG}.tar.gz"
if [ ! -d "$DIALS_INSTALLER" ]; then
  echo "Directory ${DIALS_INSTALLER} not found!"
  exit 1
fi
if [ ! -f "${DIALS_TAR}" ]; then
  tar czf ${DIALS_TAR} ${DIALS_INSTALLER}
fi
if [ ! -z "${DIALS_DEST_DIR}" ]; then
  rsync -avz ${DIALS_TAR} ${DIALS_DEST_DIR}/${DIALS_VERSION}
fi
#
# MAC GRAPHICAL INSTALLER
#
if [ "`uname`" = "Darwin" ]; then
  cd $DIALS_INSTALLER
  # slightly buried for a cleaner top-level directory, although this may cause
  # more user confusion than it prevents
  ./install --prefix=/Applications --compact --no-app
  if [ $? -ne 0 ]; then
    echo "Fatal installer error"
    exit 1
  fi
  #APP_BIN="${APP_DIR}/phenix-${DIALS_VERSION}/build/bin"
  #${APP_BIN}/libtbx.create_mac_app phenix \
  #  --app_name="Phenix-${DIALS_VERSION}" \
  #  --python_interpreter=/usr/bin/python \
  #  --dest="/Applications/Phenix-${DIALS_VERSION}" > py2app.log
  # also copy phenix_env.* files to top-level directory
  ${PYTHON_BIN} ./lib/libtbx/auto_build/create_mac_pkg.py \
    --package_name=DIALS --version=${DIALS_VERSION} \
    --background=${BUILD_DIR}/modules/dials/util/installer_background.jpg \
    --license=${BUILD_DIR}/modules/dials/license.txt \
    --organization=gov.lbl.dials /Applications/dials-${DIALS_VERSION}
  DIALS_PKG="DIALS-${DIALS_VERSION}-${HOST_TAG}.pkg.zip"
  if [ ! -z "${DIALS_DEST_DIR}" ]; then
    rsync -avz ${DIALS_PKG} ${DIALS_DEST_DIR}/${DIALS_VERSION}
  fi
  cd ..
fi
# get rid of bulky directories
#if [ $DIALS_DEBUG -eq 0 ] && [ $DEV_BUILD -eq 1 ]; then
#  rm -rf $INST_DIR
#  rm -rf ${DIALS_INSTALLER}
#fi
