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
# but some customization is recommended, e.g. via a thin wrapper script.
#
if [ -z "${DIALS_DEST_DIR}" ]; then
  # XXX This should be defined to wherever installers live, minus version
  DIALS_DEST_DIR=""
fi
if [ -z "${DIALS_DEBUG}" ]; then
  DIALS_DEBUG=0
fi
#
# Stuff below here should be site-independent...
#
DIALS_VERSION=$1
HOST_TAG=$2
if [ -z "$DIALS_VERSION" ]; then
  echo "DIALS_VERSION must be second argument!"
  exit 1
fi
if [ -z "$HOST_TAG" ]; then
  echo "HOST_TAG must be third argument!"
  exit 1
fi
BUILD_DIR="`libtbx.find_in_repositories dials`/../.."
echo "BUILD_DIR is ${BUILD_DIR}"
DEV_BUILD=`echo ${DIALS_VERSION} | grep -c dev`
PYTHON_BIN=${BUILD_DIR}/build/bin/libtbx.python
$PYTHON_BIN ${BUILD_DIR}/modules/dials/util/assemble_installer.py \
  --version=${DIALS_VERSION} --host-tag=${HOST_TAG} $BUILD_DIR
DIALS_INSTALLER="dials-installer-${DIALS_VERSION}"
DIALS_TAR="dials-installer-${DIALS_VERSION}-${HOST_TAG}.tar.gz"
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
if [ $DIALS_DEBUG -eq 0 ] && [ $DEV_BUILD -eq 1 ]; then
  rm -rf $INST_DIR
  rm -rf ${DIALS_INSTALLER}
fi
