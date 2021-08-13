##############################################
Processing Small Molecule MicroED/3DED: Biotin
##############################################

Author: Jessica Bruhn, `NanoImaging Services <https://www.nanoimagingservices.com/>`_

.. highlight:: none

.. warning::

  This tutorial was prepared using DIALS version 3.5.4 for Mac, downloaded
  from :doc:`this site <../../../installation>`. Results may differ with other
  versions of the software.

General Notes
=============

* Raw data for biotin can be downloaded from |biotin|
* Data were collected with a Glacios TEM and a CETA-D camera using
  Leginon, by the continuous rotation method.
* Information about Leginon [1]_
* Information about the sample and data collection [2]_
* Final structure deposited in CCDC:
  `BIOTIN16 <https://dx.doi.org/10.5517/ccdc.csd.cc27ydsd>`_
* Biotin: |C10H16N2O3S|, |P212121|, (5.1 Å, 10.4 Å, 20.8 Å, 90°, 90°, 90°)

.. |biotin| image:: https://zenodo.org/badge/DOI/10.3389/fmolb.2021.648603.svg
                  :target: https://zenodo.org/record/4737864#.YRYpH3VKhhE

.. [1] .. pubmed:: 33030237 Leginon

.. [2] .. pubmed:: 34327213 Small Molecule Microcrystal Electron Diffraction for the Pharmaceutical Industry

.. |C10H16N2O3S| replace:: C\ :sub:`10`\ H\ :sub:`16`\N\ :sub:`2`\O\ :sub:`3`\S

.. |P212121| replace:: P2\ :sub:`1`\ 2\ :sub:`1`\2\ :sub:`1`


Note about pedestal and offset
==============================

* Data collected with the CETA/CETA-D cameras includes some pixels for
  which the recorded intensity value is negative. Many of these pixels
  likely still contain valuable information and therefore should not be
  entirely excluded from the data.
* SMV format, a common format for X-ray crystallography diffraction images,
  does not allow for “unsigned data” (negative values). Because
  of this, Leginon adds an offset value to all pixels when converting to
  SMV format to make all values positive. Leginon also excludes negative
  values beyond a reasonable threshold.
* Some data processing programs, such as DIALS, allow for unsigned
  data. Therefore, it can be beneficial to add back in a reasonable
  pedestal value, essentially resetting the zero point.
* For data collected with Leginon, both the offset value applied
  (LEGINON_OFFSET), and a suggested pedestal value (IMAGE_PEDESTAL) are
  listed in the image metadata.

Import images
=============

Correct interpretation of these images requires a special format class,
installed as a plug-in for the dxtbx library. Installation only has to
be done once, using this command:

.. code-block:: bash

  dxtbx.install_format -u https://raw.githubusercontent.com/dials/dxtbx_ED_formats/master/FormatSMVCetaD_TUI.py

With that in place, you can start data processing. Start with dataset
801406_1.

.. code-block:: bash

    cd 801406_1
    mkdir DIALS
    cd DIALS
    dials.import ../*.img panel.pedestal=980

.. note::

    The suggested pedestal can be found by opening one of the images
    with a text editor and finding the suggested pedestal value. This
    value is provided by Leginon: |image001|.

.. |image001| image:: https://dials.github.io/images/Biotin_NIS/image001.png

Find the beam centre
====================

.. code-block:: bash

    dials.image_viewer imported.expt

* Opens the imported experiment in the Image Viewer
* Turn off |image003|
    * Note that the beam center from the image header is not accurate
* It can be helpful to find the beam center using an image with good
  diffraction spots. Try moving the slider at the top of the window to
  image 45
* Change Zoom to 50%. For weaker data you can also increase the
  brightness value
* Move the mouse to the center of the direct beam, not the center of
  the beamstop. It can be helpful to find Friedel pairs and draw lines
  between them. The beam center should be in the center of Friedel pairs.
  |image004|
* Make a note of the slow and fast beam center values at the bottom of
  the window (red box).
* Close the viewer

.. |image003| image:: https://dials.github.io/images/Biotin_NIS/image003.png

.. |image004| image:: https://dials.github.io/images/Biotin_NIS/image004.png

Re-import with the correct beam center
======================================

.. code-block:: bash

    dials.import ../*.img slow_fast_beam_centre=988,1059 panel.pedestal=980

Generate a mask for the beam center (optional)
==============================================

.. code-block:: bash

    dials.generate_mask imported.expt untrusted.circle=1059,988,100

* The first two numbers are the beam center and the second is the
  diameter of the mask. Note that the beam center values are flipped
  compared to the import step.

Find spots
==========

.. code-block:: bash

    dials.find_spots imported.expt gain=0.1 d_min=1.2 mask=pixels.mask\
      min_spot_size=2 max_separation=15 max_spot=5000

* This step can be very time consuming when working with a new
  detector/data collection parameters. You want to make sure you are
  detecting enough spots corresponding to the lattice of interest to
  index your data and should adjust parameters until this is the case.
* Here I have set the gain to 0.1, which is not the true gain for this
  detector. Modifying the gain here allows for the detection
  of more spots, but this will not impact integration step, which will
  use the correct gain (>26).
* Note that I have adjusted the min and max spot sizes, as well as
  their separation.
* I have also set the d_min to 1.2 Å to reduce the impact from high
  resolution spots not related to the lattice of interest (“zingers”).

In the log file (``dials.find_spots.log``), note the number of spots
found on each image. In this case, there are very few spots found in
the first 9 images:

.. code-block::

    Found 28 strong pixels on image 1
    Found 8 strong pixels on image 2
    Found 27 strong pixels on image 3
    Found 1 strong pixels on image 4
    Found 8 strong pixels on image 5
    Found 160 strong pixels on image 6
    Found 55 strong pixels on image 7
    Found 67 strong pixels on image 8
    Found 50 strong pixels on image 9
    Found 276 strong pixels on image 10
    Found 347 strong pixels on image 11
    Found 598 strong pixels on image 12
    Found 584 strong pixels on image 13
    Found 483 strong pixels on image 14
    Found 506 strong pixels on image 15
    Found 327 strong pixels on image 16
    Found 422 strong pixels on image 17
    Found 286 strong pixels on image 18
    Found 413 strong pixels on image 19
    Found 483 strong pixels on image 20
    Found 328 strong pixels on image 21
    Found 142 strong pixels on image 22
    Found 189 strong pixels on image 23
    Found 140 strong pixels on image 24
    Found 424 strong pixels on image 25
    Found 1031 strong pixels on image 26
    Found 735 strong pixels on image 27
    Found 419 strong pixels on image 28
    Found 374 strong pixels on image 29
    Found 500 strong pixels on image 30
    Found 670 strong pixels on image 31
    Found 718 strong pixels on image 32

This is evident when opening the images with the image viewer:

.. code-block:: bash

    dials.image_viewer imported.expt strong.refl

.. image:: https://dials.github.io/images/Biotin_NIS/image008.png

For now, just make a mental note that there are very few spots on images 1-9.

Indexing
========

.. code-block:: bash

    dials.index imported.expt strong.refl detector.fix=distance

* Fixing the detector distance is essential for electron diffraction
  data, especially given the unusual geometry relative to X-ray
  diffraction.
* Make sure that the camera length (distance) is carefully calibrated
  for your microscope as this value will not be refined by DIALS.

In the log file (``dials.index.log``), note the final ``RMSD_X`` and
``RMSD_Y``. The smaller the value the better. Generally, values lower than
3 are acceptable for electron diffraction data.

.. code-block::

    RMSDs by experiment:
    +-------+--------+----------+----------+------------+
    |   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
    |    id |        |     (px) |     (px) |   (images) |
    |-------+--------+----------+----------+------------|
    |     0 |    986 |   1.2036 |   1.8388 |    0.34601 |
    +-------+--------+----------+----------+------------+


Also note the % of spots indexed. 79% is quite good for electron
diffraction, but lower values (~30%) are still okay.

.. code-block::

    +------------+-------------+---------------+-------------+
    |   Imageset |   # indexed |   # unindexed | % indexed   |
    |------------+-------------+---------------+-------------|
    |          0 |        1101 |           299 | 78.6%       |
    +------------+-------------+---------------+-------------+

Find the Bravais lattice (optional)
===================================

.. code-block:: bash

    dials.refine_bravais_settings indexed.refl indexed.expt detector.fix_list=distance

* Potential lattices are listed.
* Note the ``Metric fit`` and ``rmsd`` values, as well as the
  recommended solutions:

  .. code-block::

    Chiral space groups corresponding to each Bravais lattice:
    aP: P1
    mP: P2 P21
    oP: P222 P2221 P21212 P212121
    +------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+----------+
    |   Solution |   Metric fit |   rmsd | min/max cc   |   #spots | lattice   | unit_cell                                |   volume | cb_op    |
    |------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+----------|
    |   *      5 |       0.4805 |  0.061 | 0.756/0.867  |     1013 | oP        | 5.53  11.00  22.27  90.00  90.00  90.00  |     1354 | a,b,c    |
    |   *      4 |       0.4805 |  0.061 | 0.793/0.793  |     1012 | mP        | 5.42  10.79  21.81  90.00  90.13  90.00  |     1277 | a,b,c    |
    |   *      3 |       0.4776 |  0.059 | 0.756/0.756  |     1002 | mP        | 5.79  23.31  11.50  90.00  90.25  90.00  |     1552 | -a,-c,-b |
    |   *      2 |       0.4495 |  0.059 | 0.867/0.867  |     1003 | mP        | 12.96   6.54  26.30  90.00  89.45  90.00 |     2230 | -b,-a,-c |
    |   *      1 |       0      |  0.062 | -/-          |      986 | aP        | 5.19  10.36  20.82  90.36  90.31  90.32  |     1120 | a,b,c    |
    +------------+--------------+--------+--------------+----------+-----------+------------------------------------------+----------+----------+
    * = recommended solution

* Lattice choice is generally less straightforward for electron
  diffraction data.
* When in doubt, process the data in P1 and determine the true symmetry
  after processing many datasets or after phasing the data (ADDSYM in
  Platon is great for doing this).
* In this case, we know that the biotin crystal should be |P212121|,
  Solution #5 (primitive orthorhombic), but let’s just process in P1 to
  start with.

Refine the cell
===============

.. code-block:: bash

    dials.refine indexed.refl indexed.expt scan_varying=True\
      detector.fix=all parameterisation.block_width=0.25\
      beam.fix="all in_spindle_plane out_spindle_plane *wavelength"\
      beam.force_static=False beam.smoother.absolute_num_intervals=1

* We fix the detector position and orientation with ``detector.fix=all``,
  but now we are allowing the crystal unit cell, crystal orientation, and
  beam direction parameters to vary on a frame-by-frame basis with
  ``scan_varying=True``.

Now the ``RMSD_X`` and ``RMSD_Y`` have decreased significantly:

.. code-block::

    RMSDs by experiment:
    +-------+--------+----------+----------+------------+
    |   Exp |   Nref |   RMSD_X |   RMSD_Y |     RMSD_Z |
    |    id |        |     (px) |     (px) |   (images) |
    |-------+--------+----------+----------+------------|
    |     0 |    903 |  0.79972 |   1.1824 |    0.23681 |
    +-------+--------+----------+----------+------------+

Integration
===========

.. code-block:: bash

    dials.integrate refined.expt refined.refl d_min=0.8

* We have set the high-resolution limit (``d_min``) to 0.8 Å. Applying a
  resolution cutoff at integration speeds up later steps, especially
  scaling multiple datasets together. You want to set this to a smaller
  d-spacing limit than you expect for your dataset.

Scaling
=======

.. code-block:: bash

    dials.scale integrated.expt integrated.refl d_min=0.8

Though we will not directly use the output from scaling individual
datasets, performing this step at this stage is helpful to assess the
quality of the individual dataset.

Find the file ``dials.scale.html`` and open it in a web browser.

* Note the useful statistics by resolution shells.
* Keep in mind that merging statistics from incomplete and low
  multiplicity data are less reliable. It is generally best to wait to
  assess the final resolution cutoff until data from multiple crystals
  have been combined.
* Scroll down the page a little and click |image019|. This brings up two
  graphs. Let’s focus on the "Scale and R\ :sub:`merge` vs batch" plot:
  |image020|
* This plots the scale factor and R\ :sub:`merge` on a per frame (N)
  basis. Let’s focus on the orange R\ :sub:`merge`  line (right axis).
* Note that there is an uptick in R\ :sub:`merge` at the beginning and
  the end of the dataset. The higher R\ :sub:`merge` values at the start
  are likely due to the low number of spots that were found on those
  images, due to suboptimal centering. The uptick at the end is more
  likely to be due to radiation damage.
* We will remove these high R\ :sub:`merge` frames after combining data
  from all four crystals.

.. |image019| image:: https://dials.github.io/images/Biotin_NIS/image019.png

.. |image020| image:: https://dials.github.io/images/Biotin_NIS/image020.png


Other datasets
==============

Repeat this process for the other three datasets.

Hint, here are the import commands I used for each dataset:

801406_1
    ``dials.import ../*.img slow_fast_beam_centre=988,1059 panel.pedestal=980``

801574_1
    ``dials.import ../*.img slow_fast_beam_centre=992,1022 panel.pedestal=831``

802003_1
    ``dials.import ../*.img slow_fast_beam_centre=986,1026 panel.pedestal=791``

810542_1
    ``dials.import ../*.img slow_fast_beam_centre=998,1024 panel.pedestal=1619``

Multi-dataset symmetry determination
====================================

We will run ``dials.cosym`` in a new directory alongside the dataset
directories

.. code-block:: bash

    mkdir cosym
    cd cosym
    dials.cosym ../801406_1/integrated.{expt,refl}\
      ../801574_1/integrated.{expt,refl}\
      ../802003_1/integrated.{expt,refl}\
      ../810542_1/integrated.{expt,refl}

Towards the end of the log we see:

.. code-block::

    Scoring all possible sub-groups
    +-------------------+----+--------------+----------+--------+--------+---------+--------------------+
    | Patterson group   |    |   Likelihood |   NetZcc |   Zcc+ |   Zcc- |   delta | Reindex operator   |
    |-------------------+----+--------------+----------+--------+--------+---------+--------------------|
    | P 1 2/m 1         | *  |        0.55  |     2.45 |   9.49 |   7.03 |     0.2 | -a,-c,-b           |
    | P m m m           |    |        0.321 |     7.93 |   7.93 |   0    |     0.3 | -a,-b,c            |
    | P -1              |    |        0.051 |    -7.93 |   0    |   7.93 |     0   | -a,-b,c            |
    | P 1 2/m 1         |    |        0.039 |    -1.31 |   7.04 |   8.35 |     0.2 | -a,-b,c            |
    | P 1 2/m 1         |    |        0.039 |    -1.33 |   7.02 |   8.35 |     0.3 | -b,-a,-c           |
    +-------------------+----+--------------+----------+--------+--------+---------+--------------------+
    Best solution: P 1 2/m 1
    Unit cell: (5.18887, 20.8461, 10.2932, 90, 90.1936, 90)
    Reindex operator: -a,-c,-b
    Laue group probability: 0.550
    Laue group confidence: 0.355
    Reindexing operators:
    x,y,z: [0, 1, 2, 3]

* Note that the suggested Patterson group is ``P 1 2/m 1``, not ``P m m m``
  as we expect for biotin.
* This is likely due to the data being fed into dials.cosym being a
  little suboptimal.

* We could force ``dials.cosym`` to choose |P212121| by adding
  ``space_group=P212121`` to the above command and move on if we
  wanted, or we could improve the input files.

Excluding bad frames
====================

Let’s improve the input files by going back and reprocess all datasets
to exclude bad frames.

* Remember the higher R\ :sub:`merge` we saw for certain frames in
  ``dials.scale.html`` for the first dataset? Let’s re-process each
  dataset and remove these “bad” frames.
* Start at the import step and exclude some of the bad frames:

801406_1
    ``dials.import ../*.img slow_fast_beam_centre=988,1059 panel.pedestal=980 image_range=6,129``

    .. image:: https://dials.github.io/images/Biotin_NIS/801406_1.png
       :width: 50%

801574_1
    ``dials.import ../*.img slow_fast_beam_centre=992,1022 panel.pedestal=831 image_range=1,101``

    .. image:: https://dials.github.io/images/Biotin_NIS/801574_1.png
       :width: 50%

802003_1
    ``dials.import ../*.img slow_fast_beam_centre=986,1026 panel.pedestal=791 image_range=1,126``

    .. image:: https://dials.github.io/images/Biotin_NIS/802003_1.png
       :width: 50%

810542_1
    ``dials.import ../*.img slow_fast_beam_centre=998,1024 panel.pedestal=1619 image_range=3,128``

    .. image:: https://dials.github.io/images/Biotin_NIS/810542_1.png
       :width: 50%

Process as before and re-run ``dials.cosym`` with the trimmed data:

.. code-block::

    +-------------------+-----+--------------+----------+--------+--------+---------+--------------------+
    | Patterson group   |     |   Likelihood |   NetZcc |   Zcc+ |   Zcc- |   delta | Reindex operator   |
    |-------------------+-----+--------------+----------+--------+--------+---------+--------------------|
    | P m m m           | *** |        0.973 |     9.49 |   9.49 |   0    |     0.5 | a,b,c              |
    | P 1 2/m 1         |     |        0.012 |     0.38 |   9.75 |   9.36 |     0.5 | -a,-c,-b           |
    | P 1 2/m 1         |     |        0.007 |    -0.18 |   9.38 |   9.55 |     0.3 | -b,-a,-c           |
    | P 1 2/m 1         |     |        0.007 |    -0.21 |   9.35 |   9.56 |     0.5 | a,b,c              |
    | P -1              |     |        0.001 |    -9.49 |   0    |   9.49 |     0   | a,b,c              |
    +-------------------+-----+--------------+----------+--------+--------+---------+--------------------+
    Best solution: P m m m
    Unit cell: (5.19177, 10.294, 20.8491, 90, 90, 90)

Now the ``P m m m`` Patterson group is the most likely, as expected.

Scale the data together
=======================

Starting from the output of ``dials.cosym``:

.. code-block:: bash

    dials.scale symmetrized.expt symmetrized.refl nproc=8 d_min=0.8

* The ``d_min=0.8`` is not actually necessary because we only
  integrated to 0.8 Å.
* Open ``dials.scale.html``

.. image:: https://dials.github.io/images/Biotin_NIS/image026.png

* Note the increase in R\ :sub:`merge` part way through collection of
  dataset #1 (801574_1).
* Let’s remove some of those images and see how that changes things:
  ``dials.scale symmetrized.expt symmetrized.refl nproc=8 d_min=0.8 exclude_images="1:49:101"``
* Here we have removed images 49-101 from dataset #1 as these had a
  fairly high R\ :sub:`merge` |image028|
* Note that the completeness in the lower resolution shells have
  decreased a small amount. Let’s try adding back in some frames to boost
  the completeness back to 100% in the low-resolution shells:
  ``dials.scale symmetrized.expt symmetrized.refl nproc=8 d_min=0.8 exclude_images="1:60:101"``
  |image030|
* This looks a lot better in terms of completeness.
* Looking at ``dials.scale.html`` we can probably improve this a little
  by removing some images from the end of dataset #2 |image032|
* So, we run ``dials.scale symmetrized.expt symmetrized.refl nproc=8 d_min=0.8 exclude_images="1:60:101" exclude_images="2:121:126"``
  |image034| |image036|
* This looks reasonably good

.. |image028| image:: https://dials.github.io/images/Biotin_NIS/image028.png

.. |image030| image:: https://dials.github.io/images/Biotin_NIS/image030.png

.. |image032| image:: https://dials.github.io/images/Biotin_NIS/image032.png

.. |image034| image:: https://dials.github.io/images/Biotin_NIS/image034.png

.. |image036| image:: https://dials.github.io/images/Biotin_NIS/image036.png

Export the data
===============

We export the scaled, unmerged dataset to MTZ format:

.. code-block:: bash

    dials.export scaled.refl scaled.expt mtz.hklout=biotin_P222.mtz

* It can be helpful to open the ``biotin_P222.mtz`` file in CCP4’s
  ViewHKL program
* You want to see a nice decay of intensities, with more intense spots
  in the middle and lower intensity spots towards the edges. You also
  want to make sure that you see some variation in your reflection
  intensities. They should not all be the same value as can happen when
  scaling goes poorly.

.. image:: https://dials.github.io/images/Biotin_NIS/image038.png

Now we want to convert to SHELX format for structure solution

.. code-block:: bash

    xia2.to_shelx biotin_P222.mtz biotin_P222 CHNOS

Solve the structure
===================

Prepare an ``.ins`` file for SHELXT or SHELXD by either running XPREP or
manually editing the .ins file as shown below:

.. code-block::

    TITL biotin_P212121 in P2(1)2(1)2(1)
    CELL 0.02501   5.19177  10.29400  20.84910  90.0000  90.0000  90.0000
    ZERR    4.00   0.00104   0.00206   0.00417   0.0000   0.0000   0.0000
    LATT -1
    SYMM 0.5-X, -Y, 0.5+Z
    SYMM -X, 0.5+Y, 0.5-Z
    SYMM 0.5+X, 0.5-Y, -Z
    SFAC C H N O S
    UNIT 40 64 8 12 4
    FIND    10
    PLOP    14    17    19
    MIND 1.0 -0.1
    NTRY 1000
    HKLF 4
    END

Now phase the data:

.. code-block:: bash

    shelxt biotin_P212121

If you are processing a more challenging organic small molecule dataset
you could try this:

.. code-block:: bash

    shelxt biotin_P212121 -y -m1000

If this fails (it should not fail for biotin as this is really
high-quality data), you could try SHELXD, which has been more
successful in phasing challenging datasets here at NIS:

.. code-block:: bash

    shelxd biotin_P212121

* Note that you have to have the space group correct for SHELXD to work.
* When you are having difficulties, try solving this in P1 and figuring
  out the proper space group once you have a solution with ADDSYM in
  PLATON.
* It can help to increase the ``NTRY``. Try 50000 for challenging cases.

If SHELXD fails, I usually go back and remove more datasets/bad images
and try again.

If that fails, you can try molecular replacement with PHASER.

* Note that you will need merged data and an R-free set. I recommend
  using ``dials.merge`` and then ``freerflag`` in CCP4.
* The model used needs to be very accurate in terms of RMSD with the
  final structure.
* When defining the contents of the ASU, try setting this to 2% solvent
  content.

If you are still struggling, I recommend going back and collecting more
data or growing better crystals. Sometimes one crystal will diffract to
much higher resolution than the others. For challenging cases, we have
collected data from ~200 crystals just to find ~5 good ones to combine.

**Good luck!**
