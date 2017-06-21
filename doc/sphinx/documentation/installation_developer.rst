+++++++++++++++++++++++++++
Installation for Developers
+++++++++++++++++++++++++++

Setting up a Development Environment on Linux or Mac
====================================================

Prerequisites:  make and change to a working directory to contain the new source code
and build. Then download these bootstrap modules::

  wget https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py

or if wget is not available::

  curl https://raw.githubusercontent.com/cctbx/cctbx_project/master/libtbx/auto_build/bootstrap.py > bootstrap.py

Then::

  python bootstrap.py --builder=dials

Explanation:  Several steps are performed: hot, update, base, build.  If desired, they can be run individually at the command line::

  python bootstrap.py --builder=dials hot
  python bootstrap.py --builder=dials update
  python bootstrap.py --builder=dials base
  python bootstrap.py --builder=dials build

* "hot" downloads static tarballs containing prerequisite dependencies, to the module directory
* "update" anonymously checks out or updates source code for dials and cctbx, to the module directory
* "base" downloads and installs python and third party python packages to the base directory
* "build" configures and compiles dials and cctbx

The DIALS builder pulls code from repositories in multiple locations. If you are a developer with an account at cci.lbl.gov you can add --cciuser=<your cci user name>. This will allow you to make commits back to the source tree. Similarly, if you have a source forge account, you can add --sfuser=<your source forge user name>. If you have an account at GitHub and have installed a registered SSH key on your machine, you can add --git-ssh.

For subsequent login sessions, be sure to set the environment in order to use the command-line dispatchers::

  source build/setpaths.sh # or setpaths.csh for tcsh

Additional packages can be installed in the modules directory, e.g., download the dials_regression tarball from the main dials web
page. Additional modules may be configured as in this example::

  libtbx.configure dials_regression

Restarting the Base Install if One Component Fails
--------------------------------------------------

It has required quite a bit of experimentation to get the "base" install correct.
Here is a procedure to restart the base install if it dies in the middle, and needs to be
restarted.  First, the top of the base output gives a list of python packages to be installed.
On linux it looks something like this::

  python numpy hdf5 biopython freetype gettext glib expat fontconfig render pixman png tiff cairo gtk fonts wxpython
  matplotlib pyopengl imaging reportlab misc

Identify the subset of packages that has failed to install; as an example assume that wxpython and subsequent packages still need to be
installed.  Then run the base installer using the just-installed python as the "with-python" base::

  python modules/cctbx_project/libtbx/auto_build/install_base_packages.py \
  --with-python=`pwd`/base/bin/python wxpython matplotlib pyopengl imaging reportlab misc

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

Building Using Your Own Custom Python
-------------------------------------

Using your own python to configure and build DIALS requires the dependencies listed above to already be installed and properly configured. Instruct the builder to use your python thusly::

  python bootstrap.py --builder=dials --with-python=<absolute path to your python binary>

This will run the 'hot' and 'update' steps as normal, skip the 'base' step that would normally download and build python and necessary DIALS python dependencies, and finally run the 'build' step with the python you specified. Note that the DIALS team doesn't recommend this procedure, but it can be necessary if using DIALS in concert with other packages.

.. _build_dials_windows:

Installing DIALS from SVN on Windows
====================================

Dependencies
------------

In order to follow this process, you will need the following programs
installed on your machine:

 - 64bit python (https://www.python.org/downloads/release/python-278/)
 - HDF5 (http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/hdf5-1.8.14-win64-vs2012-shared.zip)
 - A subversion client

Before trying to compile anything, you will need to add the path to the hdf5.h
file to the INCLUDE environment variable. If you don't have the INCLUDE
enviroment variable, just add it. You will also need a C/C++ compiler (e.g.
visual sudio).

Getting the CCTBX and DIALS sources
-----------------------------------

To get this started, create a directory to contain the cctbx build::

  mkdir cctbx
  cd cctbx

Download http://cci.lbl.gov/cctbx_build/results/current/cctbx_bundle_for_installer.tar.gz
and unpack into the directory ``cctbx\modules``.

Now checkout the cctbx sources into the ``cctbx\modules\cctbx_project`` directory::

  cd modules
  git clone https://github.com/cctbx/cctbx_project.git cctbx_project

In the ``modules`` directory of your cctbx installation, checkout the dials source
in the following way::

  git clone https://github.com/dials/dials.git dials

This may take some time, but will fetch all the dials source code and deposit in
in a folder called dials within the cctbx source directory.

Now, create a build directory in ``cctbx\build``. and configure the cctbx
installation and build the c++ libraries as follows::

  cd ..
  mkdir build
  cd build
  python ..\modules\cctbx_project\libtbx\configure.py dials
  setpaths.bat
  libtbx.scons

Note that the setpaths.bat script needs to be sourced each time you want to build
cctbx or run a cctbx program.

You should now be good to go!

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
