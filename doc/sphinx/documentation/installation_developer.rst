+++++++++++++++++++++++++++
Installation for Developers
+++++++++++++++++++++++++++

Setting up a Development Environment on Linux or Mac
====================================================

Prerequisites:  make and change to a working directory to contain the new source code
and build. Then download these bootstrap modules::

  wget https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py

or if wget is not available::

  curl https://raw.githubusercontent.com/dials/dials/master/installer/bootstrap.py > bootstrap.py

Then::

  python bootstrap.py

Explanation:  Several steps are performed: update, base, build.  If desired, they can be run individually at the command line::

  python bootstrap.py update
  python bootstrap.py base
  python bootstrap.py build

* "update" checks out or updates source code to the 'modules' directory
* "base" downloads and installs python and third party python packages to the 'conda_base' directory
* "build" configures and compiles dials and cctbx

For subsequent login sessions, be sure to set the environment in order to use the command-line dispatchers::

  source build/setpaths.sh # or setpaths.csh for tcsh

Additional packages can be installed in the modules directory, e.g., download the dials_regression tarball from the main dials web
page. Additional modules may be configured as in this example::

  libtbx.configure dials_regression

Creating a Relocatable Installer Bundle on Linux
------------------------------------------------

Starting with the developer build just created, we can create a tarball suitable for public
distribution.  Caveat is that we build our 64-bit installer on Centos 5.4, so that most conceivable users
will be installing on a more modern OS back-compatible with the installer.

Change to the working directory used above.  Then::

  ./modules/dials/installer/dials_installer.sh

..creates an installer called tmp/dials-installer-dev.tar.gz

This can be relocated to a new directory, untarred, then::

  cd dials-installer-dev
  ./install -h [prints a help message]
  ./install --prefix=[absolute path for relocated dials installation]

Downloading the DIALS regression test data
==========================================

The DIALS regression test data, needed for some of the DIALS tests, can be
obtained `here <http://dials.diamond.ac.uk/developers/dials_regression.tgz>`_::

  cd ../modules
  curl http://dials.diamond.ac.uk/developers/dials_regression.tgz > dials_regression.tgz
  tar -xzvf dials_regression.tgz
  libtbx.configure dials_regression

For those with svn access to the CCI server, it can also be obtained as
follows. Checkout the data into the cctbx source
directory and configure as follows, replacing "USERNAME" for your username::


  cd ../modules
  svn checkout svn+ssh://USERNAME@cci.lbl.gov/dials_regression/trunk dials_regression
  libtbx.configure dials_regression
