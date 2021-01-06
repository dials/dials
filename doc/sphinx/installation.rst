++++++++++++
Installation
++++++++++++

.. include:: installation.stable_release

Development Builds
==================

Nightly build installers are available for Linux and Mac OS.
For instructions on compiling from source or setting up a DIALS development
environment, see :doc:`/documentation/installation_developer`.

.. button::
   :text: Mac installer
   :link: https://dials.diamond.ac.uk/diamond_builds/dials-macosx-conda3.pkg

.. button::
   :text: Linux installer
   :link: https://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64-conda3.tar.gz


Installation
============

Mac graphical binary installers
-------------------------------

We provide a graphical package installer for Mac users. Download the
`Mac OS X <https://dials.diamond.ac.uk/diamond_builds/dials-macosx-conda3.pkg>`_
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

  wget https://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64-conda3.tar.xz
  tar -xJf dials-linux-x86_64-conda3.tar.xz
  cd dials-installer-dev

Or on Mac::

  curl https://dials.diamond.ac.uk/diamond_builds/dials-macosx-conda3.tar.gz > dials-macosx.tar.gz
  tar -xzf dials-macosx.tar.gz
  cd dials-installer-dev

Then to install in the /usr/local directory (you may need to add a ``sudo``
before the command)::

  ./install

or to install in a specified directory::

  ./install --prefix=/path/to/installation/directory/

To use DIALS, open a new terminal window and type, e.g.::

  source /path/to/installation/directory/dials-dev/dials_env.sh
