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
  echo "DIALS_VERSION must be first argument!"
  exit 1
fi
if [ -z "$HOST_TAG" ]; then
  echo "HOST_TAG must be second argument!"
  exit 1
fi
DIALS="`libtbx.find_in_repositories dials`"
MISC_OPTIONS=""
if [ ! -z "${DIALS_DEST_DIR}" ]; then
  MISC_OPTIONS="--dist-dir=${DIALS_DEST_DIR}"
fi
if [ "${DIALS_DEBUG}" != "0" ]; then
  MISC_OPTIONS="${MISC_OPTIONS} --debug"
fi
libtbx.make_dist \
  --version=${DIALS_VERSION} \
  --host-tag=${HOST_TAG} \
  ${MISC_OPTIONS} \
  ${DIALS}/util/make_dist.phil
if [ $? -ne 0 ]; then
  echo "Fatal error assembling installer"
  exit 1
fi
