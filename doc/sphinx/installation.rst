++++++++++++
Installation
++++++++++++

Stable Release: DIALS 1.1.3
===========================

The current stable releases can be downloaded from `Github <https://github.com/dials/dials/releases/tag/v1.1.0>`_,
where you can also find further `release notes <https://github.com/dials/dials/releases/tag/v1.1.0>`_.

.. button::
   :text: DIALS 1.1.3 Mac installer (OS X 10.11)
   :link: https://github.com/dials/dials/releases/download/v1.1.0/dials-v1-1-3-macosx.pkg

.. button::
   :text: DIALS 1.1.3 Mac installer (OS X 10.6)
   :link: https://github.com/dials/dials/releases/download/v1.1.0/dials-v1-1-3-macosx-10.6.pkg

.. button::
   :text: DIALS 1.1.3 Linux installer
   :link: https://github.com/dials/dials/releases/download/v1.1.0/dials-v1-1-3-linux-x86_64.tar.xz

.. button::
   :text: DIALS 1.1.3 Source installer
   :link: https://github.com/dials/dials/releases/download/v1.1.0/dials-v1-1-3-source.tar.xz

Development Builds
==================

Nightly build installers are available for Linux and Mac OS and may be
downloaded from `LBL <http://cci.lbl.gov/dials/installers/>`_ or
`Diamond <http://dials.diamond.ac.uk/diamond_builds/>`_.
Builds for Microsoft Windows are not currently available, but will be added in
the near future.
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


Windows binary installers
-------------------------

Unfortunately we don't currently provide Windows binaries, although we do plan
to add them in due course. For instructions on building DIALS from source, see
:ref:`build_dials_windows`
