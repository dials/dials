Installer
=========

Nightly build installers are available for Linux and Mac OS and may be
downloaded from http://cci.lbl.gov/dials/installers/ and
http://dials.diamond.ac.uk/builds/.  Builds for Microsoft Windows are not
currently available, but will be added in the near future.

Setting up a development environment on Linux
----------------------------------------------

Prerequisites: make and change to a working directory to contain the new source
code and build. Then download these bootstrap modules:

.. code-block:: none

  svn export svn://svn.code.sf.net/p/cctbx/code/trunk/libtbx/auto_build/bootstrap.py
  svn export svn://svn.code.sf.net/p/cctbx/code/trunk/libtbx/auto_build/bootstrap_public.py

Both codes are needed; bootstrap.py provides the basic functions;
bootstrap_public.py is a thin wrapper permitting general users to download all
the source code with anonymous access to the Berkeley lab servers.

Then:

.. code-block:: none

  python bootstrap_public.py --builder=dials

Explanation: Several steps are performed: hot, update, base, build. If desired,
they can be run individually at the command line:

.. code-block:: none

  python bootstrap_public.py --builder=dials hot
  python bootstrap_public.py --builder=dials update
  python bootstrap_public.py --builder=dials base
  python bootstrap_public.py --builder=dials build

:hot: downloads static tarballs containing prerequisite dependencies, to the module directory
:update: anonymously checks out or updates source code for dials and cctbx, to the module directory
:base: downloads and installs python and third party python packages to the base directory
:build: configures and compiles dials and cctbx

For subsequent login sessions, be sure to set the environment in order to use
the command-line dispatchers:

.. code-block:: none
  
  source build/setpaths.sh # or setpaths.csh for tcsh

Additional packages can be installed in the modules directory, e.g., download
the dials_regression tarball from the main dials web page. Additional modules
may be configured as in this example:

.. code-block:: none

  libtbx.configure dials_regression

Restarting the Base Install if One Component Fails
--------------------------------------------------

It has required quite a bit of experimentation to get the "base" install
correct. Here is a procedure to restart the base install if it dies in the
middle, and needs to be restarted. First, the top of the base output gives a
list of python packages to be installed. On linux it looks something like this:

.. code-block:: none
  
  python numpy hdf5 biopython freetype gettext glib expat fontconfig render pixman png tiff cairo gtk fonts wxpython matplotlib pyopengl imaging reportlab misc

Identify the subset of packages that has failed to install; as an example
assume that wxpython and subsequent packages still need to be installed. Then
run the base installer using the just-installed python as the "with-python"
base:

.. code-block:: none

  python modules/cctbx_project/libtbx/auto_build/install_base_packages.py \
    --with-python=`pwd`/base/bin/python \
    wxpython \
    matplotlib \
    pyopengl \
    imaging \
    reportlab \
    misc

Creating a relocatable installer bundle on Linux
------------------------------------------------

Starting with the developer build just created, we can create a tarball
suitable for public distribution. Caveat is that we build our 64-bit installer
on Centos 5.4, so that most conceivable users will be installing on a more
modern OS back-compatible with the installer.

Change to the working directory used above. Then:

.. code-block:: none

  ./modules/dials/installer/dials_installer.sh

..creates an installer called tmp/dials-installer-dev.tar.gz

This can be relocated to a new directory, untarred, then:

.. code-block:: none
  
  cd dials-installer-dev
  ./install -h [prints a help message]
  ./install --prefix=[absolute path for relocated dials installation]

Install DIALS from SVN on Linux
===============================

Dependencies
------------

In order to follow this process, you will need to following tools installed on
your machine:

   - csh
   - curl
   - perl
   - svn

You will also need C/C++ and Fortran compilers.

Installing Python
-----------------

To get this started, create a directory to contain the cctbx build and get the
auto_builder script

.. code-block:: none

   mkdir cctbx
   cd cctbx
   svn export https://svn.code.sf.net/p/cctbx/code/trunk/libtbx/auto_build

And then:

.. code-block:: none

   mkdir build
   cd build
   ../auto_build/install --all

Then let this run for some time. The following packages are installed this way:

|    biopython-1.58.tar.gz
|    freetype-2.4.2.tar.gz
|    h5py-2.0.1-edit.tar.gz
|    hdf5-1.8.8.tar.bz2
|    Imaging-1.1.7.tar.gz
|    matplotlib-1.2.0.tar.gz
|    numpy-1.6.1.tar.gz
|    PyRTF-0.45.tar.gz
|    reportlab-2.6.tar.gz
|    wxPython-src-2.8.12.1.tar.gz

though the versions of the packages may depend on the OS version you have. If
you experience problems with the auto_build script, you can either use your
system python or install python from https://www.python.org/ and install the
required packages manually.

Getting the CCTBX and DIALS sources
-----------------------------------

The cctbx provides a script to help in checking out the cctbx svn repository. To
get started with this, go back to the cctbx root directory and download
everything using this script:

.. code-block:: none

   cd ..
   svn export https://svn.code.sf.net/p/cctbx/code/trunk/libtbx/development/cctbx_svn_getting_started.csh
   ./cctbx_svn_getting_started.csh

Then wait a while again while everything downloads. You will now have a
directory named "sources" containing all the cctbx source code. In the "sources"
directory of your cctbx installation, checkout the dials source in the following
way:

.. code-block:: none

   cd sources
   svn checkout https://svn.code.sf.net/p/dials/code/trunk dials

This may take some time, but will fetch all the dials source code and deposit in
in a folder called dials within the cctbx source directory.  To include dials
within the cctbx build process, execute the following command.

.. code-block:: none

   libtbx.configure dials

Now compile the DIALS sources by executing the following commands in the "build"
directory.

.. code-block:: none

   cd ../build
   ./base/bin/python ../sources/cctbx_project/libtbx/configure.py dials
   . setpaths.sh
   make

Note that the setpaths.sh script needs to be sourced each time you want to build
dials or run a dials program; this can be added to your .bashrc file if
necessary.

You should now be good to go!

Install DIALS from SVN on Windows
=================================

Dependencies
------------

In order to follow this process, you will need to following programs installed on
your machine:

 - 64bit python (https://www.python.org/downloads/release/python-278/)
 - HDF5 (http://www.hdfgroup.org/ftp/HDF5/current/bin/windows/hdf5-1.8.14-win64-vs2012-shared.zip)
 - A subversion client

Before trying to compile anything, you will need to add the path to the hdf5.h
file to the INCLUDE environment variable. If you don't have the INCLUDE
enviroment variable, just add it. You will also need a C/C++ compiler (e.g.
visual sudio).

Getting the CCTBX and DIALS sources
-----------------------------------

To get this started, create a directory to contain the cctbx build.

.. code-block:: none

   mkdir cctbx
   cd cctbx

Download http://cci.lbl.gov/cctbx_build/results/current/cctbx_bundle_for_installer.tar.gz
and unpack into the directory "cctbx\sources".

Now checkout the cctbx sources into the "cctbx\sources\cctbx_project" directory.

.. code-block:: none
   cd sources
   svn checkout svn://svn.code.sf.net/p/cctbx/code/trunk cctbx_project

In the "sources" directory of your cctbx installation, checkout the dials source
in the following way:

.. code-block:: none

   svn checkout https://svn.code.sf.net/p/dials/code/trunk dials

This may take some time, but will fetch all the dials source code and deposit in
in a folder called dials within the cctbx source directory.
  
Now, create a build directory in "cctbx\build". and configure the cctbx
installation and build the c++ libraries as follows.

.. code-block:: none

   cd ..
   mkdir build
   cd build
   python ..\sources\cctbx_project\libtbx\configure.py dials
   setpaths.bat
   libtbx.scons

Note that the setpaths.bat script needs to be sourced each time you want to build
cctbx or run a cctbx program.

You should now be good to go!

Downloading the DIALS regression test data
==========================================

To obtain the dials regression test data, needed for some of the dials tests,
you will need access to the CCI server. Checkout the data into the cctbx source
directory and configure as follows, replacing "USERNAME" for your username:

.. code-block:: none

   cd ../sources
   svn checkout svn+ssh://USERNAME@cci.lbl.gov/dials_regression/trunk dials_regression
   libtbx.configure dials_regression
