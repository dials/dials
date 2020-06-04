Switch the DIALS development environment to use conda compilers

Newly created development environments will use native conda-forge compilers rather than system compilers, and enable the use of C++11.

A new 'dials' script, located above the build directory, replaces the existing 'setpaths'-family of scripts.
This script loads the conda environment, and makes all commands from within that environment available.

