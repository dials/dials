++++++++++++
Installation
++++++++++++

.. include:: installation.stable_release

Development Builds
==================

Nightly build installers are available for Linux and Mac OS and may be
downloaded from `LBL <http://cci.lbl.gov/dials/installers/>`_ or
`Diamond <http://dials.diamond.ac.uk/diamond_builds/>`_.
Builds for Microsoft Windows are experimental and may not work as expected.
For instructions on compiling from source or setting up a DIALS development
environment, see :doc:`/documentation/installation_developer`.

.. button::
   :text: Mac installer (OS X 10.11)
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-macosx.pkg

.. button::
   :text: Mac installer (OS X 10.6)
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-macosx-10.6.pkg

.. button::
   :text: Linux installer
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64.tar.xz

.. button::
   :text: Windows archive
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-windows.zip

.. button::
   :text: Source installer
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-source.tar.xz

Installation
============

Mac graphical binary installers
-------------------------------

We provide a graphical package installer for Mac users. Download the
`Mac OS X 10.11 <http://dials.diamond.ac.uk/diamond_builds/dials-macosx.pkg>`_
or
`Mac OS X 10.6 <http://dials.diamond.ac.uk/diamond_builds/dials-macosx-10.6.pkg>`_
installer and double click the ``.pkg`` file to start the
graphical installer. Follow the instructions, which will install DIALS in the
``/Applications/`` directory. To use DIALS, open a new terminal window and type,
e.g.::

  source /Applications/dials-dev20150903/dials_env.sh

to setup the DIALS environment.


Mac and Linux binary installers
-------------------------------

We provide binary ``tar.gz`` and ``tar.xz`` files for various Mac and Linux
platforms, e.g. on Linux::

  wget http://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64.tar.xz
  tar -xJf dials-linux-x86_64.tar.xz
  cd dials-installer-dev

Or on Mac::

  curl http://dials.diamond.ac.uk/diamond_builds/dials-macosx.tar.gz > dials-macosx.tar.gz
  tar -xzf dials-macosx.tar.gz
  cd dials-installer-dev

Then to install in the /usr/local directory (you may need to add a ``sudo``
before the command)::

  ./install

or to install in a specified directory::

  ./install --prefix=/path/to/installation/directory/

To use DIALS, open a new terminal window and type, e.g.::

  source /path/to/installation/directory/dials-dev/dials_env.sh


Windows installation
--------------------

DIALS support on Windows is currently experimental and we do not provide Windows binaries.
We plan to support Windows in the near future and add binary installers in due course.

To use DIALS you need to unpack the .zip archive, open a command prompt,
navigate to the dials-installer directory, and run::

  dials_env

For instructions on building DIALS from source, see
:ref:`build_dials_windows`
