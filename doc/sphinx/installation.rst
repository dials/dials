++++++++++++
Installation
++++++++++++

Binary Installers
=================

Nightly build installers are available for Linux and Mac OS and may be
downloaded from `LBL <http://cci.lbl.gov/dials/installers/>`_ or
`Diamond <http://dials.diamond.ac.uk/diamond_builds/>`_.
Builds for Microsoft Windows are not currently available, but will be added in
the near future.
For instructions on compiling from source or setting up a DIALS development
environment, see :doc:`/documentation/installation_developer`.

Download the latest nightly build:
----------------------------------

.. button::
   :text: Mac installer (OS X 10.10.3)
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-macosx-10.10.3.tar.gz

.. button::
   :text: Mac installer (OS X 10.6)
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-macosx-10.6.tar.gz

.. button::
   :text: Linux installer
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64.tar.gz

.. button::
   :text: Source installer
   :link: http://dials.diamond.ac.uk/diamond_builds/dials-source.tar.gz


Mac and Linux binary installers
-------------------------------

We provide binary ``tar.gz`` files for various Mac and Linux platforms, e.g.::

  curl http://dials.diamond.ac.uk/diamond_builds/dials-linux-x86_64.tar.gz > dials-linux-x86_64.tar.gz
  tar -xzf dials-linux-x86_64.tar.gz
  cd dials-installer-dev

Then to install in the /usr/local directory (you may need to add a ``sudo``
before the command)::

  ./install

or to install in a specified directory::

  ./install --prefix=/path/to/installation/directory/

To use DIALS, open a new terminal window and type, e.g.::

  source /path/to/installation/directory/dials-dev/dials_env.sh


Mac graphical binary installers
-------------------------------

We also  provide a graphical package installer for Mac users. Download, e.g., the
file ``dials-dev-239-mac-intel-osx-x86_64.pkg.zip``, double click the ``.pkg.zip``
file in Finder to unzip, then once again double click the ``.pkg`` to start the
graphical installer. Follow the instructions, which will install DIALS in
``/Applications/`` directory. To use DIALS, open a new terminal window and type,
e.g. ::

  source /Applications/dials-dev-239/dials_env.sh

to setup the DIALS environment.


Windows binary installers
-------------------------

Unfortunately we don't currently provide Windows binaries, although we do plan
to add them in due course. For instructions on building DIALS from source, see
:ref:`build_dials_windows`
