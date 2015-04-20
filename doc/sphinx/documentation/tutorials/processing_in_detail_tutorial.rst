Processing in Detail
====================

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, exporting to MTZ) or you can run the
whole lot through :doc:`dials.process </programs/dials_process>`, which just
chains them together (and incidentally does all of the processing in P1). In
this tutorial we will run through each of the steps in turn, checking the output
as we go. We will also enforce the correct lattice symmetry.

Tutorial data
-------------

The following example uses a Thaumatin dataset collected using beamline I04
at Diamond Light Source which is available for download from |thaumatin|.

.. |thaumatin| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10271.png
               :target: http://dx.doi.org/10.5281/zenodo.10271

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that the image headers are read, and a file describing
their contents (:ref:`datablock.json <datablock-json>`) is written. It's worth noting that if
this file is changed subsequent processing (even with :doc:`dials.process </programs/dials_process>`) can
use this.

::

  dials.import data/th_8_2_0*cbf

The output just describes what the software understands of the images it was
passed, in this case one sweep of data containing 540 images.

::

  The following parameters have been modified:

  input {
    datablock = <image files>
  }

  --------------------------------------------------------------------------------
  DataBlock 0
    format: <class 'dxtbx.format.FormatCBFMiniPilatusDLS6MSN100.FormatCBFMiniPilatusDLS6MSN100'>
    num images: 540
    num sweeps: 1
    num stills: 0
  --------------------------------------------------------------------------------
  Writing datablocks to datablock.json

Find Spots
^^^^^^^^^^

The first "real" task in any DIALS processing will be the spot finding.
Here we request multiple processors to speed up the spot-finding (:samp:`nproc=4`).
It takes a little while because we are finding spots on every image in the
dataset. This reflects the modular philosophy of the DIALS toolkit and will
enable us to do global refinement later on.

::

  dials.find_spots datablock.json nproc=4

This will just report the number of spots found.

::

  The following parameters have been modified:

  spotfinder {
    mp {
      nproc = 4
    }
  }
  input {
    datablock = datablock.json
  }


  Configuring spot finder from input parameters
  --------------------------------------------------------------------------------
  Finding strong spots in imageset 0
  --------------------------------------------------------------------------------

  Finding spots in image 0 to 540...
  Extracting strong pixels from images (may take a while)
  Extracted strong pixels from images
  Merging 4 pixel lists
  Merged 4 pixel lists with 922120 pixels
  Extracting spots
  Extracted 219125 spots
  Calculating 219125 spot centroids
  Calculated 219125 spot centroids
  Calculating 219125 spot intensities
  Calculated 219125 spot intensities
  Filtering 219125 spots by number of pixels
  Filtered 116321 spots by number of pixels
  Filtering 116321 spots by peak-centroid distance
  Filtered 116082 spots by peak-centroid distance

  --------------------------------------------------------------------------------
  Saving 116082 reflections to strong.pickle
  Saved 116082 reflections to strong.pickle
  Time Taken: 28.706979

The default parameters for :doc:`dials.find_spots </programs/dials_find_spots>`
usually do a good job
for Pilatus images, such as these. However they may not be optimal for data from
other detector types, such as CCDs or image plates. Issues with incorrectly
set gain or sigma thresholds might lead to far too many spots being extracted
(for example). If you are having issues with spot finding, it is worth
inspecting the images with :program:`dials.image_viewer`::

  dials.image_viewer datablock.json

Viewing the various images from 'image' to 'threshold' gives an idea of how the
various parameters affect the spot finding algorithm. The final image,
'threshold' is the one on which spots are found, so ensuring this produces peaks
at real diffraction spot positions will give the best chance of success.

Having found strong spots it is worth checking the image viewer again::

  dials.image_viewer datablock.json strong.pickle

The :program:`dials.image_viewer` tool is not as fast as tools such as ADXV,
however it does integrate well with DIALS data files. Information about
the beam centre, spot centroids, reflection shoeboxes and other data stored in
the pickle files created by DIALS programs can be overlayed on the diffraction
images. You may need to adjust the colour scheme and brightness to get the best
out of it. A brightness of 20 with the 'invert' colour scheme works well with
this data. Move forward a few images to find a spot whose complete rocking curve
is recorded. The highest valued pixel in that three dimensional spot is marked
with a pink dot. The spot centre of mass is a red cross. This is usually close to
the peak pixel, but slightly offset as the centroid algorithm allows to calculate
the spot centre at a better precision than the pixel size and image angular 'width'.
The strong pixels marked as being part of the peak are highlighted with a green
dot. The reflection shoebox you see with a blue border is the smallest
three dimensional box that
can contain the continuous peak region, that is, there is no background border
region displayed here.

.. image:: /figures/found_spot.png

Indexing
^^^^^^^^

The next step will be indexing of the strong spots, by default using a 3D FFT
algorithm, although the 1D FFT algorithm can be selected using the parameter
:samp:`indexing.method=fft1d`. We will pass in all the strong spots found in
the dataset - so no need to select subsets of images widely separated in
:math:`\phi`.

::

  dials.index datablock.json strong.pickle

If known, the space group and unit cell can be
provided at this stage using the :samp:`space_group` and :samp:`unit_cell`
parameters, otherwise indexing and refinement will be carried out in the
primitive lattice using space group P1.

::

  The following parameters have been modified:

  input {
    datablock = datablock.json
    reflections = strong.pickle
  }

  Found max_cell: 199.1 Angstrom
  Setting d_min: 3.88778695204
  FFT gridding: (256,256,256)
  Number of centroids used: 13282
  model 1 (13227 reflections):
  Crystal:
      Unit cell: (57.604, 58.025, 149.904, 89.691, 89.825, 89.768)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2615, -0.9013},
                  { 0.8912,  0.3922,  0.2278},
                  { 0.2939, -0.8819,  0.3685}}
      B matrix:  {{ 0.0174,  0.0000,  0.0000},
                  {-0.0001,  0.0172,  0.0000},
                  {-0.0001, -0.0001,  0.0067}}
      A = UB:    {{ 0.0061, -0.0044, -0.0060},
                  { 0.0154,  0.0067,  0.0015},
                  { 0.0051, -0.0152,  0.0025}}


  41 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 1)
  ################################################################################


  Summary statistics for observations matched to predictions:
  -------------------------------------------------------------------------
  |                   | Min     | Q1      | Med      | Q3        | Max    |
  -------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5997 | -0.3171 | -0.2353  | -0.1522   | 0.6011 |
  | Yc - Yo (mm)      | -1.184  | -0.2289 | -0.09807 | -0.008457 | 0.3844 |
  | Phic - Phio (deg) | -0.7938 | -0.1812 | -0.0362  | 0.06264   | 0.7968 |
  | X weights         | 253.5   | 401.7   | 404.2    | 405.1     | 405.6  |
  | Y weights         | 260.5   | 402.6   | 404.5    | 405.3     | 405.6  |
  | Phi weights       | 416.1   | 527.9   | 531      | 532.2     | 533.3  |
  -------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.26289  | 0.17719  | 0.1821   |
  | 1    | 4049 | 0.042615 | 0.045835 | 0.096025 |
  | 2    | 4049 | 0.038331 | 0.043512 | 0.053594 |
  | 3    | 4049 | 0.036176 | 0.042359 | 0.032431 |
  | 4    | 4049 | 0.034808 | 0.040253 | 0.023344 |
  | 5    | 4049 | 0.033692 | 0.037341 | 0.018083 |
  | 6    | 4049 | 0.032377 | 0.034582 | 0.016841 |
  | 7    | 4049 | 0.031501 | 0.033056 | 0.016667 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.18314 | 0.19219 | 0.11112  |
  ---------------------------------------------
  Increasing resolution to 2.9 Angstrom
  model 1 (31986 reflections):
  Crystal:
      Unit cell: (57.805, 57.782, 150.024, 90.003, 89.981, 89.990)
      Space group: P 1
      U matrix:  {{ 0.3454, -0.2591, -0.9020},
                  { 0.8915,  0.3909,  0.2291},
                  { 0.2933, -0.8832,  0.3660}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000,  0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}


  121 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 2)
  ################################################################################


  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.2821 | -0.02348 | -0.002085  | 0.03725 | 0.2655 |
  | Yc - Yo (mm)      | -0.7087 | -0.02207 | -0.001341  | 0.01982 | 0.2813 |
  | Phic - Phio (deg) | -1.053  | -0.01024 | -0.0003781 | 0.01023 | 0.9065 |
  | X weights         | 229.3   | 399.6    | 403.3      | 404.8   | 405.6  |
  | Y weights         | 210.9   | 399.5    | 403.3      | 404.8   | 405.6  |
  | Phi weights       | 399.5   | 526.6    | 530.5      | 532.1   | 533.3  |
  --------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.04046  | 0.035062 | 0.016443 |
  | 1    | 4049 | 0.039277 | 0.034359 | 0.016479 |
  | 2    | 4049 | 0.039221 | 0.034184 | 0.016431 |
  | 3    | 4049 | 0.039096 | 0.033931 | 0.016391 |
  | 4    | 4049 | 0.038792 | 0.033558 | 0.016366 |
  | 5    | 4049 | 0.038204 | 0.032931 | 0.016353 |
  | 6    | 4049 | 0.037521 | 0.03224  | 0.01632  |
  | 7    | 4049 | 0.037146 | 0.03194  | 0.016279 |
  | 8    | 4049 | 0.037058 | 0.031931 | 0.01627  |
  | 9    | 4049 | 0.037049 | 0.031938 | 0.016269 |
  | 10   | 4049 | 0.037049 | 0.031939 | 0.016269 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  --------------------------------------------
  | Exp | Nref | RMSD_X | RMSD_Y  | RMSD_Z   |
  |     |      | (px)   | (px)    | (images) |
  --------------------------------------------
  | 0   | 4049 | 0.2154 | 0.18569 | 0.10846  |
  --------------------------------------------
  Increasing resolution to 1.9 Angstrom
  model 1 (91948 reflections):
  Crystal:
      Unit cell: (57.815, 57.779, 150.032, 90.018, 89.995, 89.987)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2589, -0.9020},
                  { 0.8914,  0.3909,  0.2293},
                  { 0.2932, -0.8833,  0.3658}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000,  0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}


  309 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 3)
  ################################################################################


  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.4729 | -0.02289 | 0.008137   | 0.0387  | 0.6873 |
  | Yc - Yo (mm)      | -1.421  | -0.02119 | 0.002495   | 0.02561 | 1.281  |
  | Phic - Phio (deg) | -1.434  | -0.01349 | -0.0007786 | 0.01238 | 0.9047 |
  | X weights         | 179.9   | 383.7    | 397        | 403.2   | 405.6  |
  | Y weights         | 171     | 377.3    | 394.4      | 402.5   | 405.6  |
  | Phi weights       | 318.9   | 520.1    | 529.5      | 533.3   | 533.3  |
  --------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.048284 | 0.043291 | 0.024544 |
  | 1    | 4049 | 0.046364 | 0.042412 | 0.024789 |
  | 2    | 4049 | 0.046319 | 0.042245 | 0.024662 |
  | 3    | 4049 | 0.046228 | 0.041887 | 0.024539 |
  | 4    | 4049 | 0.04606  | 0.041289 | 0.024444 |
  | 5    | 4049 | 0.045965 | 0.040565 | 0.024374 |
  | 6    | 4049 | 0.04602  | 0.040007 | 0.024329 |
  | 7    | 4049 | 0.046084 | 0.039749 | 0.024309 |
  | 8    | 4049 | 0.046109 | 0.039694 | 0.024306 |
  | 9    | 4049 | 0.046114 | 0.039689 | 0.024305 |
  | 10   | 4049 | 0.046114 | 0.039688 | 0.024305 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  --------------------------------------------
  | Exp | Nref | RMSD_X | RMSD_Y  | RMSD_Z   |
  |     |      | (px)   | (px)    | (images) |
  --------------------------------------------
  | 0   | 4049 | 0.2681 | 0.23075 | 0.16204  |
  --------------------------------------------
  Increasing resolution to 0.9 Angstrom
  model 1 (114690 reflections):
  Crystal:
      Unit cell: (57.814, 57.785, 150.036, 90.014, 89.991, 89.988)
      Space group: P 1
      U matrix:  {{ 0.3454, -0.2589, -0.9020},
                  { 0.8914,  0.3909,  0.2292},
                  { 0.2933, -0.8833,  0.3658}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000,  0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}


  342 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 4)
  ################################################################################


  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5349 | -0.03342 | -0.004698  | 0.03183 | 0.6688 |
  | Yc - Yo (mm)      | -1.409  | -0.02779 | -0.0004348 | 0.02336 | 1.265  |
  | Phic - Phio (deg) | -1.424  | -0.01409 | -6.861e-05 | 0.0144  | 0.9046 |
  | X weights         | 135.2   | 371.6    | 393.4      | 402.6   | 405.6  |
  | Y weights         | 153.1   | 361.8    | 389        | 401.2   | 405.6  |
  | Phi weights       | 318.9   | 519.5    | 530.4      | 533.3   | 533.3  |
  --------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.049701 | 0.047766 | 0.026101 |
  | 1    | 4049 | 0.049518 | 0.047234 | 0.025985 |
  | 2    | 4049 | 0.04946  | 0.046987 | 0.025939 |
  | 3    | 4049 | 0.049369 | 0.046646 | 0.025902 |
  | 4    | 4049 | 0.049268 | 0.046296 | 0.025895 |
  | 5    | 4049 | 0.049185 | 0.046    | 0.025901 |
  | 6    | 4049 | 0.049136 | 0.045853 | 0.025917 |
  | 7    | 4049 | 0.049117 | 0.045808 | 0.025936 |
  | 8    | 4049 | 0.049113 | 0.045798 | 0.025944 |
  | 9    | 4049 | 0.049113 | 0.045797 | 0.025945 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.28554 | 0.26626 | 0.17297  |
  ---------------------------------------------
  Final refined crystal models:
  model 1 (114690 reflections):
  Crystal:
      Unit cell: (57.802, 57.779, 150.018, 90.010, 89.993, 89.989)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2590, -0.9020},
                  { 0.8914,  0.3909,  0.2292},
                  { 0.2932, -0.8832,  0.3659}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000,  0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}

  Saving refined experiments to experiments.json
  Saving refined reflections to indexed.pickle


It is worth looking through this output to understand what the indexing program
has done. Note that this log (minus the preamble about modified parameters)
is automatically captured in the file :file:`dials.index.log`. There is also
a somewhat more information written into :file:`dials.index.debug.log`, but
this is probably only helpful if something has gone wrong and you are trying
to track down why.

Inspecting the log shows that the indexing step is done at fairly low
resolution: ``Setting d_min: 3.88778695204``. The resolution limit of data that
can be used in indexing is determined by the size of the 3D FFT grid and the
likely maximum cell dimension. Here we
used :math:`256^3` grid points: ``FFT gridding: (256,256,256)``.
What follows are four macrocycles
of refinement at increasing resolution to bootstrap the indexing solution to as
many of the strong reflections as possible. In each case you can see that only
4049 reflections are used in the refinement job. The diffraction geometry is
here described by only 16 parameters (6 for the detector, 1 beam angle, 3
crystal 'misset' angles and 6 triclinic cell parameters). The problem is thus
hugely overdetermined. In order to save time, refinement uses a subset of the
input reflections, by default using 50 reflections for every degree of the scan.

Continuing to look through the log, we see that the first macrocyle of refinement makes
a big improvement, reducing the positional RMSDs from 0.26 to 0.03 mm in X and
0.17 to 0.03 mm in Y. The second macrocycle includes more reflections, after
extending to 2.9 Angstroms. The current model now shows slightly worse RMSDs
at the start, now that the higher resolution reflections are included, but refinement reduces
these from 0.040 to 0.037 mm in X and from 0.035 to 0.032 mm in Y.
A similar situation is observed when the resolution is extended again to 1.9 Angstroms.
The RMSDs start higher again, now that more reflections are included, but refinement
is able to drive these down a little, from 0.048 to 0.046 mm in X and 0.043 to 0.040 mm in Y.
The final macrocycle includes data out to 0.9 Angstroms, which is well beyond
the highest resolution recorded 'strong' spot. Nevertheless, more reflections
have been included so refinement acts again to improve the models slightly
including this higher-resolution information. The final model provides
rmsds of 0.049 mm in X, 0.046 mm in Y and 0.026 degrees in :math:`\phi`, corresponding
to 0.29 pixels in X, 0.27 pixels in Y and 0.17 image widths in :math:`\phi`.

Despite the high quality of this data, we notice from the ``Summary statistics``
tables that there there are some outliers appearing as resolution increases.
In the final macrocyle we see the
distribution of positional residuals in the Y direction is tight around the
median, except for extreme values both positive and negative of more than 1 mm.
The angular residuals show a similar pattern with half the data having residuals
of less than about 0.14 degrees from the predicted positions, but the extreme
is as much as 1.4 degrees from the predicted diffraction angle. We are happy
with the indexing solution though and will deal with these outliers in the
separate refinement step to come later.

If you want to specify the Bravais lattice for processing (i.e. include the
lattice constraints in the refinement) then you need to either specify this
lattice at this stage as

::

  space_group=P4

as a command-line option to :doc:`dials.index </programs/dials_index>`
or you can use
:doc:`dials.refine_bravais_settings </programs/dials_refine_bravais_settings>`,
which will take the results of the P1 autoindexing and run refinement with all
of the possible Bravais settings applied - after which you may select the
preferred solution.

::

  dials.refine_bravais_settings experiments.json indexed.pickle

gives a table containing the metric fit, rmsds (in mm) and unit cell for
each Bravais setting...

::

  The following parameters have been modified:

  input {
    experiments = experiments.json
    reflections = indexed.pickle
  }

  -------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd #spots  crystal_system                                 unit_cell  volume      cb_op
  -------------------------------------------------------------------------------------------------------------
         9  0.0258 dg 0.070   4049   tetragonal tP  57.78  57.78 149.99  90.00  90.00  90.00  500719      a,b,c
         8  0.0258 dg 0.069   4049 orthorhombic oC  81.72  81.73 150.01  90.00  90.00  90.00 1001925  a-b,a+b,c
         7  0.0148 dg 0.068   4049 orthorhombic oP  57.78  57.76 149.98  90.00  90.00  90.00  500612      a,b,c
         6  0.0231 dg 0.068   4049   monoclinic mC  81.73  81.74 150.03  90.00  89.99  90.00 1002292  a-b,a+b,c
         5  0.0258 dg 0.069   4049   monoclinic mC  81.73  81.72 150.01  90.00  89.99  90.00 1001933 a+b,-a+b,c
         4  0.0131 dg 0.068   4049   monoclinic mP  57.77  57.79 150.00  90.00  90.01  90.00  500757   -b,-a,-c
         3  0.0148 dg 0.068   4049   monoclinic mP  57.79  57.77 150.00  90.00  89.99  90.00  500789      a,b,c
         2  0.0120 dg 0.068   4049   monoclinic mP  57.77 150.00  57.79  90.00  89.99  90.00  500842      b,c,a
         1  0.0000 dg 0.067   4049    triclinic aP  57.80  57.78 150.02  90.01  89.99  89.99  501027      a,b,c
  -------------------------------------------------------------------------------------------------------------
  usr+sys time: 0.87 seconds
  wall clock time: 8.55 seconds

In this example we would continue processing (i.e. proceed to the refinement
step, perhaps) with :samp:`bravais_setting_9.json`. Sometimes it may be
necessary to reindex the :ref:`indexed.pickle <reflection_pickle>` file output by dials.index.
However, in this case as the change of basis operator to the chosen setting
is the identity operator (:samp:`a,b,c`) this step is not needed. We run it
anyway to demonstrate its use::

  dials.reindex indexed.pickle change_of_basis_op=a,b,c

This outputs the file :ref:`reindexed_reflections.pickle <reflection_pickle>` which should be
used as input to downstream programs in place of :ref:`indexed.pickle <reflection_pickle>`.


Refinement
^^^^^^^^^^

Although the model is already refined during indexing we can also add an
explicit refinement
step using :doc:`dials.refine </programs/dials_refine>` in here. This
dataset is of exceptional quality and we are keen to squeeze the best possible
results from it. During indexing we saw the presence of outliers that we would
like to exclude from refinement, and we also used a subset of reflections. Now
we will repeat using all indexed reflections in the dataset and with outlier
rejection switched on.

As an aside, to show all the options up to and including ``expert_level = 1``
use this command::

  dials.refine -c -e 1

Equivalent command-line options exist for all the main DIALS programs.
Now, our refinement job is specified as::

  dials.refine bravais_setting_9.json reindexed_reflections.pickle \
  outlier.algorithm=tukey use_all_reflections=true

The main product of this are the files ``refined_experiments.json`` and
``refined.pickle``.

::

  The following parameters have been modified:

  refinement {
    reflections {
      use_all_reflections = true
      outlier {
        algorithm = null *tukey
      }
    }
  }
  input {
    experiments = bravais_setting_9.json
    reflections = reindexed_reflections.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  -------------------------------------------------------------------------
  |                   | Min     | Q1       | Med       | Q3      | Max    |
  -------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5281 | -0.0342  | -0.002737 | 0.03274 | 0.6464 |
  | Yc - Yo (mm)      | -1.393  | -0.02842 | 0.0005028 | 0.03071 | 1.246  |
  | Phic - Phio (deg) | -1.415  | -0.01382 | 0.0009439 | 0.0159  | 0.9129 |
  | X weights         | 135.2   | 371.6    | 393.4     | 402.6   | 405.6  |
  | Y weights         | 153.1   | 361.8    | 389       | 401.2   | 405.6  |
  | Phi weights       | 318.9   | 519.5    | 530.4     | 533.3   | 533.3  |
  -------------------------------------------------------------------------

  7002 reflections have been flagged as outliers

  Summary statistics for observations matched to predictions:
  ---------------------------------------------------------------------------
  |                   | Min      | Q1       | Med       | Q3      | Max     |
  ---------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1346  | -0.03315 | -0.002133 | 0.03256 | 0.1331  |
  | Yc - Yo (mm)      | -0.1171  | -0.02688 | 0.00085   | 0.02992 | 0.1194  |
  | Phic - Phio (deg) | -0.05839 | -0.01266 | 0.00123   | 0.01547 | 0.06047 |
  | X weights         | 135.2    | 375.2    | 394.4     | 402.7   | 405.6   |
  | Y weights         | 153.1    | 366.4    | 390.6     | 401.5   | 405.6   |
  | Phi weights       | 318.9    | 519.4    | 530.1     | 533.3   | 533.3   |
  ---------------------------------------------------------------------------

  Performing refinement...

  Refinement steps:
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 106245 | 0.046664 | 0.042696 | 0.021981 |
  | 1    | 106245 | 0.046834 | 0.042375 | 0.022009 |
  | 2    | 106245 | 0.046871 | 0.042322 | 0.021985 |
  | 3    | 106245 | 0.046909 | 0.042277 | 0.021947 |
  | 4    | 106245 | 0.046939 | 0.042241 | 0.02191  |
  | 5    | 106245 | 0.046971 | 0.042198 | 0.021891 |
  | 6    | 106245 | 0.047    | 0.042163 | 0.021885 |
  | 7    | 106245 | 0.047011 | 0.042152 | 0.021884 |
  | 8    | 106245 | 0.047012 | 0.04215  | 0.021883 |
  --------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 106245 | 0.27333 | 0.24506 | 0.14589  |
  -----------------------------------------------
  Saving refined experiments to refined_experiments.json
  Updating predictions for indexed reflections
  Saving reflections with updated predictions to refined.pickle

The effectiveness of outlier rejection can be seen from the second summary
statistics table. Now the positional residuals are all within 0.14 mm and the
worst angular residual is just 0.06 degrees. After removing reflections too
close to the spindle and doing outlier rejection, refinement still has
106245 reflections to work with, amounting to 93% of the reflections in
:file:`reindexed_reflections.pickle`.

We have done the best we can with a static model for the experiment. However,
a better model for the crystal might allow small misset rotations to occur
over the course of the scan. There are usually even small changes to the
cell dimensions (typically resulting in a net increase in cell volume) caused
by exposure to radiation during data collection. To account for both of these
effects we can extend our parameterisation to obtain a smoothed 'scan-varying'
model for both the crystal orientation and unit cell. To do this, we run a
further refinement job starting from the output of the previous job::

  dials.refine refined_experiments.json refined.pickle \
  outlier.algorithm=tukey use_all_reflections=true  \
  scan_varying=true output.experiments=sv_refined_experiments.json

Note we also overrode the default experiments output filename to avoid
overwriting the output of the earlier scan-static job. Refinement output for
this job is::

  The following parameters have been modified:

  output {
    experiments = sv_refined_experiments.json
  }
  refinement {
    parameterisation {
      crystal {
        scan_varying = true
      }
    }
    reflections {
      use_all_reflections = true
      outlier {
        algorithm = null *tukey
      }
    }
  }
  input {
    experiments = refined_experiments.json
    reflections = refined.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5273 | -0.03485 | -0.003623  | 0.03184 | 0.6418 |
  | Yc - Yo (mm)      | -1.4    | -0.02935 | -0.0006238 | 0.02901 | 1.246  |
  | Phic - Phio (deg) | -1.41   | -0.01482 | -0.0002073 | 0.01475 | 0.9106 |
  | X weights         | 135.2   | 371.6    | 393.4      | 402.6   | 405.6  |
  | Y weights         | 153.1   | 361.8    | 389        | 401.2   | 405.6  |
  | Phi weights       | 318.9   | 519.5    | 530.4      | 533.3   | 533.3  |
  --------------------------------------------------------------------------

  7215 reflections have been flagged as outliers

  Summary statistics for observations matched to predictions:
  ---------------------------------------------------------------------------
  |                   | Min      | Q1       | Med        | Q3      | Max    |
  ---------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1348  | -0.03374 | -0.003056  | 0.03165 | 0.1318 |
  | Yc - Yo (mm)      | -0.1168  | -0.02759 | -0.0001738 | 0.02832 | 0.1165 |
  | Phic - Phio (deg) | -0.05917 | -0.01363 | 3.539e-05  | 0.01426 | 0.0591 |
  | X weights         | 135.2    | 375.3    | 394.5      | 402.8   | 405.6  |
  | Y weights         | 153.1    | 366.6    | 390.6      | 401.5   | 405.6  |
  | Phi weights       | 318.9    | 519.4    | 530.1      | 533.3   | 533.3  |
  ---------------------------------------------------------------------------

  Performing refinement...

  Refinement steps:
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 106028 | 0.046722 | 0.042096 | 0.021847 |
  | 1    | 106028 | 0.046434 | 0.039706 | 0.021777 |
  | 2    | 106028 | 0.046445 | 0.039653 | 0.021761 |
  | 3    | 106028 | 0.046452 | 0.039618 | 0.021707 |
  | 4    | 106028 | 0.04646  | 0.039596 | 0.02158  |
  | 5    | 106028 | 0.046465 | 0.039588 | 0.021416 |
  | 6    | 106028 | 0.046466 | 0.039587 | 0.021335 |
  | 7    | 106028 | 0.046466 | 0.039587 | 0.021326 |
  | 8    | 106028 | 0.046466 | 0.039587 | 0.021325 |
  --------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 106028 | 0.27015 | 0.23016 | 0.14217  |
  -----------------------------------------------
  Saving refined experiments to sv_refined_experiments.json
  Updating predictions for indexed reflections
  Saving reflections with updated predictions to refined.pickle

In this case we didn't alter the default choices that affect scan-varying
refinement, the most important of which is the number of intervals into which
the full scan is divided. This determines the number of samples that will be
used by the Gaussian smoother. More samples allows sharper changes to the model,
but overdoing this will lead to unphysical changes to the model that are just
fitting noise in the data. Figuring out the optimum number of points to use
is challenging. Here we are happy with the default interval width of 36 degrees
(this is a parameter at ``expert_level = 1``).

To view the smoothly varying crystal cell parameters use the following command::

  dials.plot_scan_varying_crystal sv_refined_experiments.json

This program creates a directory :file:`scan-varying_crystal` containing
plots :file:`orientation.png` and :file:`unit_cell.png`. The latter of these
is useful to check that changes to the cell during processing appear reasonable.

.. image:: /figures/unit_cell.png

We see an overall increase in all three cell parameters, however the greatest
change, in lengths *a* and *b*, is only about 0.02 Angstroms. If
significant cell volume increases had been observed that might be indicative of
radiation damage. However we can't yet conclude that there is *no* radiation
damage from the *lack* of considerable change observed. We can at least see from
this and the low final refined RMSDs that this is a very well-behaved dataset
though.

Integration
^^^^^^^^^^^

After the refinement is done the next step is integration, which is performed
by the program :doc:`dials.integrate </programs/dials_integrate>`. Mostly, the
default parameters are fine, which will perform XDS-like 3D profile fitting. However,
for datasets with very weak background, such as this, the default :samp:`nsigma`
background outlier rejection algorithm tends to underestimate the real background
value. This is because that method is only really appropriate for values from
a normal distribution, which is a poor approximation for a Poisson distibution
with a small mean, and significant skewness. For this reason we switch off
all outlier rejection from the background calculation.

From checking the output of :samp:`dials.integrate -c` we see that the full
parameter to do this is given by :samp:`integration.background.simple.outlier.algorithm=null`
but partial string matching can be used for command line parameters when the
partial match is unambiguous. This saves a lot of typing!

We will also increase the number of processors used to speed the job up.

::

  dials.integrate sv_refined_experiments.json refined.pickle \
  outlier.algorithm=null nproc=4


Checking the log output we see that after loading in the reference reflections
from :file:`refined.pickle`,
new predictions are made up to the highest resolution at the corner of the
detector. This is fine, but if we wanted to we could have adjusted the
resolution limits using parameters :samp:`dmin` and :samp:`dmax`. The predictions
are made using the scan-varying crystal model recorded in
:file:`sv_refined_experiments.json`. This ensures that prediction is made using
the smoothly varying lattice and orientation that we determined in the refinement
step. As this scan-varying model was determined in advance of integration, each
of the integration jobs is independent and we can take advantage of true
parallelism during processing.

The profile model is then calculated from the reflections in
:file:`refined.pickle`. First reflections with a too small 'zeta'
factor are filtered out. This essentially removes reflections that are too
close to the spindle axis. In general these reflections require significant
Lorentz corrections and as a result have less trustworthy intensities anyway.
From the remaining reflection shoeboxes, the average beam divergence and
reflecting range is calculated, providing the two Guassian width parameters
:math:`\sigma_D` and :math:`\sigma_M` used in the 3D profile model.

Following this, the independent integration jobs are set up. These jobs overlap,
so reflections are assigned to one or more jobs. What follows are blocks of
information specific to each integration job.

After these jobs are finished, the reflections are 'post-processed', which includes
the application of the LP correction to the intensities. Then summary tables
are printed giving quality statistics first by frame, and then by resolution bin.
The latter of these tables and the final overall summary are reproduced here::

  Summary vs resolution
  ------------------------------------------------------------------------------------------------------
  ID | d min | # full | # part | # over | # ice | # sum | # prf | <Ibg> | <I/sigI> | <I/sigI> | <CC prf>
     |       |        |        |        |       |       |       |       |  (sum)   |  (prf)   |
  ------------------------------------------------------------------------------------------------------
  0  | 1     | 351    | 3      | 0      | 0     | 263   | 199   | 0.04  | 0.37     | 0.53     | 0.09
  0  | 1     | 1122   | 5      | 0      | 0     | 996   | 872   | 0.04  | 0.43     | 0.53     | 0.08
  0  | 1     | 2479   | 13     | 0      | 0     | 2215  | 2018  | 0.05  | 0.52     | 0.59     | 0.09
  0  | 1     | 4061   | 27     | 0      | 0     | 3658  | 3474  | 0.05  | 0.56     | 0.69     | 0.11
  0  | 1     | 5882   | 33     | 0      | 0     | 5276  | 5060  | 0.05  | 0.59     | 0.76     | 0.12
  0  | 1     | 8000   | 46     | 0      | 0     | 7100  | 6850  | 0.06  | 0.65     | 0.84     | 0.14
  0  | 1     | 10549  | 61     | 0      | 0     | 9364  | 9112  | 0.06  | 0.78     | 0.98     | 0.17
  0  | 1     | 13780  | 82     | 0      | 0     | 12256 | 11951 | 0.07  | 0.91     | 1.14     | 0.20
  0  | 1     | 18748  | 107    | 0      | 0     | 16727 | 16373 | 0.07  | 0.99     | 1.23     | 0.22
  0  | 1     | 22396  | 162    | 0      | 0     | 20130 | 19919 | 0.08  | 1.20     | 1.48     | 0.25
  0  | 1     | 26177  | 694    | 0      | 0     | 23687 | 23373 | 0.09  | 1.46     | 1.76     | 0.29
  0  | 1     | 27776  | 885    | 0      | 0     | 24489 | 24354 | 0.09  | 1.75     | 2.09     | 0.34
  0  | 1     | 28056  | 820    | 0      | 0     | 25781 | 25667 | 0.10  | 2.16     | 2.53     | 0.39
  0  | 1     | 28135  | 785    | 0      | 0     | 24629 | 24518 | 0.12  | 2.67     | 3.08     | 0.45
  0  | 1     | 28267  | 755    | 0      | 0     | 25090 | 24975 | 0.14  | 3.48     | 3.92     | 0.52
  0  | 1     | 28388  | 859    | 0      | 0     | 25955 | 25834 | 0.18  | 4.80     | 5.29     | 0.59
  0  | 1     | 28456  | 932    | 0      | 0     | 24994 | 24889 | 0.24  | 6.51     | 7.05     | 0.66
  0  | 2     | 28745  | 813    | 0      | 0     | 26045 | 25942 | 0.28  | 8.76     | 9.26     | 0.69
  0  | 2     | 28814  | 803    | 0      | 0     | 25523 | 25419 | 0.34  | 12.62    | 13.07    | 0.71
  0  | 3     | 28982  | 1096   | 0      | 0     | 26204 | 26081 | 0.41  | 25.20    | 25.26    | 0.72
  ------------------------------------------------------------------------------------------------------

  Summary for experiment 19
  ----------------------------------------------------------------
  Item                                  | Overall | Low    | High
  ----------------------------------------------------------------
  dmin                                  | 1.17    | 1.47   | 1.17
  dmax                                  | 151.25  | 151.25 | 1.47
  number fully recorded                 | 369164  | 281796 | 87368
  number partially recorded             | 8981    | 8442   | 539
  number with overloaded pixels         | 0       | 0      | 0
  number in powder rings                | 0       | 0      | 0
  number processed with summation       | 330382  | 252397 | 77985
  number processed with profile fitting | 326880  | 251052 | 75828
  <ibg>                                 | 0.17    | 0.20   | 0.07
  <i/sigi> (summation)                  | 5.61    | 7.07   | 0.91
  <i/sigi> (profile fitting)            | 6.00    | 7.46   | 1.13
  <cc prf>                              | 0.46    | 0.54   | 0.19
  cc_pearson sum/prf                    | 1.00    | 1.00   | 0.89
  cc_spearman sum/prf                   | 0.96    | 0.98   | 0.77
  ----------------------------------------------------------------

Graphical analysis of the output
^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^

Much more information is available from the integration output in graphical form
using the command

::

  dials.analyse_output integrated.pickle

By default the plots will be written into a new directory :file:`analysis` with
subdirectories for different types of analysis::

  analysis
  ├── background
  ├── centroid
  ├── intensity
  ├── reference
  └── strong

Some of the most useful plots are

* :file:`background/background_model_mean_vs_xy.png`, which shows the mean
  background value as a function of detector position.

* :file:`centroid/centroid_mean_diff_vs_phi.png`, which shows how the average
  residuals in each of X, Y, and :math:`\phi` vary as a fuction of :math:`\phi`.
  If scan-varying refinement has been successful in capturing the real changes
  during the scan then we would expect these plots to be straight lines.

  .. image:: /figures/centroid_mean_diff_vs_phi.png

* :file:`centroid/centroid_xy_residuals.png`, on which the X, Y residuals are shown
  directly. The key point here is to look for a globular shape centred at 0.0.

  .. image:: /figures/centroid_xy_residuals.png

* :file:`centroid/centroid_diff_x.png` and :file:`centroid/centroid_diff_y.png`,
  which show the difference between predicted and observed reflection positions
  in either X or Y as functions of detector position. From these plots it is very
  easy to see whole tiles that are worse than their neighbours, and either whether
  those tiles might be simply shifted or slightly rotated compared to the model
  detector.

  .. image:: /figures/centroid_diff_x.png

  .. image:: /figures/centroid_diff_y.png

* :file:`reference/reflection_corr_vs_xy.png` and
  :file:`reference/reference_corr_vs_xy.png`. These are useful companions to the
  plots of centroid residual as a function of detector position displayed above.
  Whereas the earlier plots show systematic errors in the positions and
  orientations of tiles of a multi-panel detector, these plots indicate what
  effect that (and any other position-specific systematic error) has on the
  integrated data quality. The first of these plots shows the correlation
  between reflections and their reference profiles for all reflections in the
  dataset. The second shows only the correlations between the strong reference
  reflections and their profiles (thus these are expected to be higher and do
  not extend to such high resolution). The first plot is probably the most
  useful, and that is reproduced here.

  .. image:: /figures/reflection_corr_vs_xy.png

* :file:`intensity/ioversigma_vs_z.png`. This reproduces the
  :math:`\frac{I}{\sigma_I}` information versus frame number given in the log
  file in a graphical form. Here we see that :math:`\frac{I}{\sigma_I}` is fairly
  flat over the whole dataset, which we might use as an indication that there
  were no bad frames, not much radiation damage occurred and that scale factors
  are likely to be fairly uniform.

  .. image:: /figures/ioversigma_vs_z.png


Exporting as MTZ
^^^^^^^^^^^^^^^^

The final step of dials processing is to export the integrated results to mtz
format, suitable for input to downstream processing programs such as pointless_
and aimless_.

::

  dials.export_mtz integrated.pickle refined_experiments.json hklout=integrated.mtz

And this is the output, showing the reflection file statistics.

::

  The following parameters have been modified:

  hklout = integrated.mtz
  input {
    experiments = refined_experiments.json
    reflections = integrated.pickle
  }

  Removing 23974 reflections with negative variance
  Removing 27291 profile reflections with negative variance
  Removing 2 reflections with I/Sig(I) < -5.0
  Removing 0 profile reflections with I/Sig(I) < -5.0
  Removing 4034 incomplete reflections
  Title: from dials.export_mtz
  Space group symbol from file: P4
  Space group number from file: 75
  Space group from matrices: P 4 (No. 75)
  Point group symbol from file: 4
  Number of batches: 540
  Number of crystals: 1
  Number of Miller indices: 322844
  Resolution range: 149.997 1.16993
  History:
  Crystal 1:
    Name: XTAL
    Project: DIALS
    Id: 1
    Unit cell: (57.7818, 57.7818, 149.997, 90, 90, 90)
    Number of datasets: 1
    Dataset 1:
      Name: FROMDIALS
      Id: 1
      Wavelength: 0.97625
      Number of columns: 14
      label        #valid  %valid    min     max type
      H            322844 100.00%   0.00   46.00 H: index h,k,l
      K            322844 100.00%   0.00   47.00 H: index h,k,l
      L            322844 100.00%   0.00  114.00 H: index h,k,l
      M_ISYM       322844 100.00%   1.00    8.00 Y: M/ISYM, packed partial/reject flag and symmetry number
      BATCH        322844 100.00%   2.00  539.00 B: BATCH number
      IPR          322844 100.00%  -1.79 2860.10 J: intensity
      SIGIPR       322844 100.00%   0.04   53.51 Q: standard deviation
      I            322844 100.00% -24.60 3059.53 J: intensity
      SIGI         322844 100.00%   0.09   55.45 Q: standard deviation
      FRACTIONCALC 322844 100.00%   1.00    1.00 R: real
      XDET         322844 100.00%   6.54 2456.33 R: real
      YDET         322844 100.00%   5.78 2520.61 R: real
      ROT          322844 100.00%  82.01  162.69 R: real
      LP           322844 100.00%   0.00    0.76 R: real


What to do Next
---------------

The following demonstrates how to take the output of dials processing and
continue with downstream analysis using pointless_ to sort the data and assign
the correct symmetry, followed by scaling with aimless_ and intensity analysis
using ctruncate_::

  pointless hklin integrated.mtz hklout sorted.mtz > pointless.log
  aimless hklin sorted.mtz hklout scaled.mtz > aimless.log << eof
  resolution 1.3
  anomalous off
  eof
  ctruncate -hklin scaled.mtz -hklout truncated.mtz \
  -colin '/*/*/[IMEAN,SIGIMEAN]' > ctruncate.log

to get merged data for downstream analysis. The output from this includes
the merging statistics which will give a better idea about data quality. It is
easiest to view these logfiles using the program :program:`logview`, e.g.::

  logview aimless.log

Often passing in a sensible resolution limit to aimless is helpful. Here we
assumed we ran first without a resolution limit to help decide where to cut
the data. This indicated slightly anisotropic diffraction, with diffraction along
the *c*\* direction a little better than *a*\* and *b*\* directions, which are
equivalent. Diffraction quality is good, however completeness falls off sharply,
especially in the *c*\* direction. Following this we chose to exclude all data
at a resolution higher than 1.3 Angstroms, to ensure about 80% completeness in
the outer shell. Here is the summary from aimless.log:

::

  Summary data for        Project: DIALS Crystal: XTAL Dataset: FROMDIALS

                                             Overall  InnerShell  OuterShell
  Low resolution limit                      150.00    150.00      1.32
  High resolution limit                       1.30      6.88      1.30

  Rmerge  (within I+/I-)                     0.063     0.024     0.421
  Rmerge  (all I+ and I-)                    0.071     0.027     0.489
  Rmeas (within I+/I-)                       0.077     0.030     0.580
  Rmeas (all I+ & I-)                        0.079     0.030     0.610
  Rpim (within I+/I-)                        0.044     0.017     0.398
  Rpim (all I+ & I-)                         0.034     0.014     0.357
  Rmerge in top intensity bin                0.029        -         -
  Total number of observations              307606      2538      5907
  Total number unique                        62350       552      2647
  Mean((I)/sd(I))                             10.7      27.0       1.5
  Mn(I) half-set correlation CC(1/2)         0.999     0.999     0.708
  Completeness                                98.2      99.8      80.6
  Multiplicity                                 4.9       4.6       2.2

  Anomalous completeness                      92.3     100.0      48.1
  Anomalous multiplicity                       2.4       3.0       1.5
  DelAnom correlation between half-sets      0.013     0.187    -0.003
  Mid-Slope of Anom Normal Probability       0.954       -         -

.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html
