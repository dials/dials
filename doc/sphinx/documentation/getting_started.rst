+++++++++++++++
Getting started
+++++++++++++++

Installation
============

You can download the latest release of DIALS from :doc:`this site
<../../../installation>` for Linux or Mac. Unless you intend to develop
with the DIALS toolkit, we recommend downloading the latest Stable
Release.

DIALS can be built on Windows, but this is somewhat experimental. See
`the wiki <https://github.com/dials/dials/wiki/DIALS-on-Windows>`_ for
details.

The installers are prepared as compressed archives. Once the download
has completed, you should unpack using the appropriate commands on your
system. For example, on Linux:

.. code-block:: bash

    tar xJf dials-v<MAJOR>-<MINOR>-<PATCH>-linux-x86_64.tar.xz

The values of ``<MAJOR>``, ``<MINOR>`` and ``<PATCH>`` will have to be
substituted appropriately for the version of DIALS you have downloaded.

Once this completes, a new directory, ``dials-installer`` will be
created. The installer script requires that it is run within that
directory. For a "system-wide" installation in ``/usr/local``, you
can run this in an environment where you have the appropriate write
permissions:

.. code-block:: bash

    cd dials-installer
    ./install

Alternatively you can install for the current user in a ``dials/``
directory under the home directory as follows:

.. code-block:: bash

    cd dials-installer
    ./install --prefix=$HOME/dials


Sourcing the DIALS environment
==============================

Once installation is complete, to use the software we must put the
suite of DIALS programs in the system ``PATH``. This is achieved by
sourcing a script. For example, in the BASH shell, and assuming DIALS
was installed in ``$HOME/dials``:

.. code-block:: bash

    source $HOME/dials/dials-v<MAJOR>-<MINOR>-<PATCH>/dials_env.sh

As before, the values of ``<MAJOR>``, ``<MINOR>`` and ``<PATCH>`` will
have to be substituted appropriately for the version of DIALS you have
installed.

There is also a ``dials_env.csh`` script for those who use a C shell
rather than a Bourne-like shell.

.. note::
    Sourcing this script does not produce any output; it simply sets
    the environment appropriately for DIALS programs.

To check that DIALS is correctly set up, you may like to try:

.. code-block:: bash

    dials.version

If everything has worked properly you should see some output giving
a DIALS version number, a Python version, and your installation path.

Built-in help
=============

The major DIALS command line programs typically print a useful help
message if run without any arguments. For example, try

.. code-block:: bash

    dials.import

This provides examples of usage, a list of command-line options and a
help message describing the function of the program. All DIALS programs
support the option

.. code-block:: bash

    dials.program -c

which will print a structured description of the parameters that the
program will accept. By default, only the basic parameters are shown.
To display all parameters up to an ``expert_level`` of 2, you would
enter

.. code-block:: bash

    dials.program -c -e2

In addition, it can be useful to display an expected type, a help
string and other useful information about the parameters. We do that by
increasing the ``attributes_level``.

.. code-block:: bash

    dials.program -c -e2 -a2

Parameters
==========

Apart from the command-line switches, all DIALS programs also accept
parameters in the form ``parameter=value``. In most cases this will be
sufficient though some less frequently used options may require "name
space" clarification e.g. `index_assignment.method=local`. More complex
parameter specifications can be written into a file, say
``myparams.phil`` and passed into the DIALS program as an input file.

Next steps
==========

At this point you are ready to start processing data! We recommend
checking out the :doc:`tutorials <tutorials/index>` for further details.

