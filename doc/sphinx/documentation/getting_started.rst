Getting started
===============

The aim of this tutorial is to provide a quick start guide to running DIALS for
integrating good quality synchrotron X-ray diffraction data. Many caveats apply
to this:

 - your mileage may vary etc. as the program is in development
 - we don't promise that the data are as good as from e.g. XDS
 - DIALS only does data processing you will need to use e.g. Aimless for the
   subsequent scaling

That said this tutorial illustrates that the program can work and give sensible
results.

Introduction
------------

The philosophy behind DIALS is to be explicit in performing the various steps of
data analysis rather than giving one big tool (though one big tool is available)
- DIALS is first and foremost a toolkit for doing the data analysis to be used
within other systems.

dials.process
-------------

In the simplest case, :doc:`dials.process </programs/dials_process>`
``/here/are/all/images*.cbf`` will do sensible processing, with a static model
of the experiment and sample, and will output a reflection file integrated.mtz
containing the intensity measurements assuming everything works correctly.
Some sensible options to explore are:

 - :samp:`scan_varying=true` - allow the crystal orientation and unit cell
   constants to vary during the scan
 - :samp:`mp.nproc=1` - only use one processor [#f1]_
 - :samp:`profile.fitting=False` - force summation integration instead of
   the default profile fitting
 - :samp:`block_size=N` - for some N, split the data set into N degree blocks
   for integration, so as not to overload the computer. A sensible default will
   be chosen, but use this to override that choice.
 - :samp:`-i` - pass the images to process through the standard input e.g. from
   :samp:`find . -name *.cbf` to avoid issues with limited command-line lengths

The workflow (if you do things step by step,) is as
follows:

 - import the data - this reads the image headers & makes sense of the data,
   forming sweeps etc.
 - find spots - as it says on the tin, can adjust spot size
 - index - index these strong spots, can use a variety of methods
 - refine - refine the experimental geometry (included in indexing for
   scan-static refinement)
 - integrate - actually integrate the data
 - export - create MTZ files

Some detail of the options for each of these will follow below. How they are
used can be seen a complete example script, found
:download:`here<../user-tutorial/tutorial.sh>`, which can be run as follows::

  ./tutorial.sh /path/to/data

Import
------

There are two ways of importing images

.. code-block:: none

  dials.import /path/to/the/data/*.img

or

.. code-block:: none

  find /path/to/the/data -name *.img | dials.import -i

The latter of these is useful when you have a limited number of command-line
options. This creates datablock.json which is a DIALS description of this data
set.

Find Spots
----------

Most useful parameter here is the minimum spots size. By default it is 6, but
this can be overridden with min_spot_size=N where N is e.g. 3. This takes
datablock.json and creates strong.pickle:

.. code-block:: none

  dials.find_spots min_spot_size=3 datablock.json

(for example)

Index
-----

The indexing for DIALS offers a substantial number of options - these are
detailed in the Phil file for dials.index, which is shown when you run the
program. The most useful options are:

.. code-block:: none

  unit_cell=a,b,c,al,be,ga
  space_group=P4 (say)
  indexing.method=fft3d (say)

The indexing works as

.. code-block:: none

  dials.index strong.pickle datablock.json [options]

and creates an "experiment" file experiment.json which details the crystal
lattices found and indexed.pickle, which is a copy of strong.pickle with Miller
indices and reflection predictions added.

If you are unsure of the symmetry and would like to know how different lattices
look in the refinement, run dials.index not specifying the symmetry and then run
dials.refine_bravais_settings.

Refinement
----------

The indexing includes refinement - if you do not wish to use a time varying
crystal model you can go straight to integration. If you do want to use a time
varying model, you will need to rerun the refinement with this new model as

.. code-block:: none

  dials.refine experiments.json indexed.pickle scan_varying=true

which will generate refined_experiments.json - this you pass on to integration.

Integration
-----------

As may be expected the integration in DIALS offers the greatest range of user
options, to control how the background is determined (including outlier pixels
in the background determination) the reflection profile parameters (used to
define the reflection mask, and by default discovered automatically) and the
actual algorithm to be used for peak integration e.g. sum3d or fft3d.

.. code-block:: none

  dials.integrate outlier.algorithm=null refined_experiments.json indexed.pickle

This reads the indexed reflections to determine strong reflections for profile
fitting and integrates the data in refined_experiments.json, using the default
background determination with no outlier rejection and XDS-style 3D profile
fitting. These commands are most likely to change and can be viewed by running

Export
------

If you have got this far everything else is easy: export the data as MTZ then
run pointless and aimless to resort and scale the data viz:

.. code-block:: none

  dials.export_mtz integrated.pickle refined_experiments.json
  pointless hklin integrated.mtz hklout sorted.mtz
  aimless hklin sorted.mtz hklout scaled.mtz

For details on pointless and aimless please refer to the CCP4 documentation.


.. rubric:: Footnotes

.. [#f1] Currently necessary for data in HDF5 files
