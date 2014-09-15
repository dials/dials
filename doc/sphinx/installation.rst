Installation
============

We are in the process of refining the build process (Nat Echols @ LBL is doing
the work for this) such that we can start with a clean directory and with a few
scripts bootstrap a full working configuration. 

This essentially needs to be done in two phases:

   - build your base python
   - build CCTBX on top

Dependencies
------------

In order to follow this process, you will need to following tools installed on
your machine:

   - csh
   - curl
   - perl
   - svn

You will also need C/C++ and Fortran compilers.

Building your Python
--------------------

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

Then let this run for some time... not sure for the moment if --all is the best
option - it does install a bunch of stuff - but it should be a good place to
start. The following packages are installed this way:

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

though the versions of the packages may depend on the OS version you have.
Certainly OS X.7 gives a different package list to the venerable RHEL5. Once
this is complete the next step is to check out all of the CCTBX components on
top...

Building CCTBX with this Python
-------------------------------

To get started with this, go back to the cctbx root directory and download
everything using this script:

.. code-block:: none

   cd ..
   svn export https://svn.code.sf.net/p/cctbx/code/trunk/libtbx/development/cctbx_svn_getting_started.csh
   ./cctbx_svn_getting_started.csh

Then wait a while again while everything downloads. You will now have a
directory named "sources" containing all the cctbx source code. To build cctbx,
go to the build directory build it.

.. code-block:: none

   cd build
   ./base/bin/python ../sources/cctbx_project/libtbx/configure.py cctbx rstbx iotbx
   . setpaths.sh
   make

This is for the basic cctbx - no labelit etc. which needs to be added and built.

Note that the setpaths.sh script needs to be sourced each time you want to build
cctbx or run a cctbx program; this can be added to your .bashrc file if
necessary.

Building DIALS with this CCTBX
------------------------------

In the "sources" directory of your cctbx installation, checkout the dials source
in the following way:

.. code-block:: none

   cd ../sources/
   svn checkout https://svn.code.sf.net/p/dials/code/trunk dials

This may take some time, but will fetch all the dials source code and deposit in
in a folder called dials within the cctbx source directory.

To include dials within the cctbx build process, execute the following command.

.. code-block:: none

   libtbx.configure dials

Then navigate to the cctbx build directory and build the dials source code.

.. code-block:: none

   cd ../build
   make

You should now be good to go!

Obtaining the DIALS regression test data
----------------------------------------

To obtain the dials regression test data, needed for some of the dials tests,
you will need access to the CCI server. Checkout the data into the cctbx source
directory and configure as follows, replacing "USERNAME" for your username:

.. code-block:: none

   cd ../sources
   svn checkout svn+ssh://USERNAME@cci.lbl.gov/dials_regression/trunk dials_regression
   libtbx.configure dials_regression
