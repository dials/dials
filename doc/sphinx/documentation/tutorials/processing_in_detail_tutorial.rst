Processing in Detail
====================

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, exporting to MTZ) or you can run
:samp:`xia2 -dials`, which makes informed choices for you at each stage. In
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
this file is changed subsequent processing can
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

  Setting spotfinder.filter.min_spot_size=3
  Configuring spot finder from input parameters
  --------------------------------------------------------------------------------
  Finding strong spots in imageset 0
  --------------------------------------------------------------------------------

  Finding spots in image 1 to 540...
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
  Found 1 possible hot spots
  Found 1 possible hot pixel(s)
  Filtering 219125 spots by number of pixels
  Filtered 116321 spots by number of pixels
  Filtering 116321 spots by peak-centroid distance
  Filtered 116082 spots by peak-centroid distance

  --------------------------------------------------------------------------------
  Saving 116082 reflections to strong.pickle
  Saved 116082 reflections to strong.pickle
  Time Taken: 71.761816

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
the pickle files created by DIALS programs can be overlaid on the diffraction
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

The next step will be indexing of the strong spots, which by default uses a 3D FFT
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

  Found max_cell: 199.0 Angstrom
  Setting d_min: 3.89
  FFT gridding: (256,256,256)
  Number of centroids used: 13298
  model 1 (13245 reflections):
  Crystal:
      Unit cell: (57.927, 58.377, 150.126, 89.734, 89.796, 89.680)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2607, -0.9015},
                  { 0.8911,  0.3923,  0.2281},
                  { 0.2942, -0.8821,  0.3679}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0001,  0.0171,  0.0000},
                  {-0.0001, -0.0001,  0.0067}}
      A = UB:    {{ 0.0060, -0.0044, -0.0060},
                  { 0.0153,  0.0067,  0.0015},
                  { 0.0051, -0.0151,  0.0025}}


  43 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 1)
  ################################################################################


  Summary statistics for 12585 observations matched to predictions:
  ------------------------------------------------------------------------
  |                   | Min     | Q1      | Med      | Q3       | Max    |
  ------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.6282 | -0.3151 | -0.2315  | -0.1485  | 0.5033 |
  | Yc - Yo (mm)      | -0.8444 | -0.2148 | -0.1216  | -0.03768 | 0.3894 |
  | Phic - Phio (deg) | -0.8988 | -0.2314 | -0.06123 | 0.04194  | 1.035  |
  | X weights         | 253.5   | 401.7   | 404.2    | 405.1    | 405.6  |
  | Y weights         | 260.5   | 402.6   | 404.5    | 405.3    | 405.6  |
  | Phi weights       | 416.1   | 527.9   | 531      | 532.2    | 533.3  |
  ------------------------------------------------------------------------

  2448 reflections have been flagged as outliers

  Summary statistics for 10137 observations matched to predictions:
  ------------------------------------------------------------------------
  |                   | Min     | Q1      | Med      | Q3       | Max    |
  ------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5793 | -0.307  | -0.2343  | -0.1609  | 0.1994 |
  | Yc - Yo (mm)      | -0.477  | -0.2013 | -0.1174  | -0.02053 | 0.3894 |
  | Phic - Phio (deg) | -0.404  | -0.1352 | -0.01158 | 0.05198  | 0.3607 |
  | X weights         | 272.3   | 401.8   | 404.2    | 405.1    | 405.6  |
  | Y weights         | 260.5   | 402.6   | 404.5    | 405.3    | 405.6  |
  | Phi weights       | 416.1   | 527.7   | 530.8    | 532.1    | 533.3  |
  ------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.26351  | 0.17489  | 0.12937  |
  | 1    | 8099 | 0.04359  | 0.048009 | 0.066931 |
  | 2    | 8099 | 0.038391 | 0.043164 | 0.04162  |
  | 3    | 8099 | 0.03623  | 0.041487 | 0.027362 |
  | 4    | 8099 | 0.034786 | 0.038974 | 0.019906 |
  | 5    | 8099 | 0.033455 | 0.036678 | 0.015428 |
  | 6    | 8099 | 0.031924 | 0.034589 | 0.014346 |
  | 7    | 8099 | 0.031019 | 0.03339  | 0.014189 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 8099 | 0.18035 | 0.19413 | 0.094595 |
  ---------------------------------------------
  Using d_min_step 0.7
  Increasing resolution to 3.2 Angstrom
  model 1 (23310 reflections):
  Crystal:
      Unit cell: (57.811, 57.789, 150.043, 89.996, 89.973, 89.988)
      Space group: P 1
      U matrix:  {{ 0.3453, -0.2591, -0.9020},
                  { 0.8915,  0.3909,  0.2290},
                  { 0.2933, -0.8832,  0.3660}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000, -0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}


  121 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 2)
  ################################################################################


  Summary statistics for 22493 observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.2816 | -0.02505 | -0.005411  | 0.02864 | 0.2515 |
  | Yc - Yo (mm)      | -0.7145 | -0.02142 | -0.0004552 | 0.01996 | 0.2835 |
  | Phic - Phio (deg) | -1.051  | -0.01046 | -0.0006884 | 0.0096  | 0.905  |
  | X weights         | 243.1   | 400.9    | 403.8      | 404.9   | 405.6  |
  | Y weights         | 239.1   | 401.3    | 404        | 405     | 405.6  |
  | Phi weights       | 401.5   | 527.6    | 530.8      | 532.2   | 533.3  |
  --------------------------------------------------------------------------

  1653 reflections have been flagged as outliers

  Summary statistics for 20840 observations matched to predictions:
  -----------------------------------------------------------------------------
  |                   | Min      | Q1       | Med        | Q3       | Max     |
  -----------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1169  | -0.0223  | -0.003617  | 0.03081  | 0.1198  |
  | Yc - Yo (mm)      | -0.09175 | -0.02192 | -0.001916  | 0.01659  | 0.08541 |
  | Phic - Phio (deg) | -0.04497 | -0.01023 | -0.0006471 | 0.009448 | 0.04195 |
  | X weights         | 243.1    | 401.1    | 403.8      | 404.9    | 405.6   |
  | Y weights         | 260.5    | 401.4    | 404        | 405      | 405.6   |
  | Phi weights       | 419      | 527.6    | 530.8      | 532.2    | 533.3   |
  -----------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.037562 | 0.028801 | 0.014065 |
  | 1    | 8099 | 0.03615  | 0.028519 | 0.014034 |
  | 2    | 8099 | 0.036101 | 0.02844  | 0.013976 |
  | 3    | 8099 | 0.03601  | 0.02835  | 0.013919 |
  | 4    | 8099 | 0.035814 | 0.028251 | 0.013892 |
  | 5    | 8099 | 0.035406 | 0.02811  | 0.013878 |
  | 6    | 8099 | 0.034861 | 0.028049 | 0.013852 |
  | 7    | 8099 | 0.034523 | 0.028151 | 0.013828 |
  | 8    | 8099 | 0.034439 | 0.028218 | 0.01382  |
  | 9    | 8099 | 0.034431 | 0.028227 | 0.013819 |
  | 10   | 8099 | 0.034431 | 0.028228 | 0.013818 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 8099 | 0.20018 | 0.16411 | 0.092123 |
  ---------------------------------------------
  Increasing resolution to 2.6 Angstrom
  model 1 (43346 reflections):
  Crystal:
      Unit cell: (57.803, 57.774, 150.032, 90.013, 89.986, 89.985)
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


  137 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 3)
  ################################################################################


  Summary statistics for 42351 observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5012 | -0.02516 | 0.006709   | 0.0326  | 0.2609 |
  | Yc - Yo (mm)      | -0.704  | -0.01313 | 0.006151   | 0.02701 | 0.2916 |
  | Phic - Phio (deg) | -1.052  | -0.01075 | -0.0004882 | 0.01031 | 0.9074 |
  | X weights         | 202.8   | 396.9    | 402.3      | 404.5   | 405.6  |
  | Y weights         | 210.9   | 396.1    | 402.2      | 404.5   | 405.6  |
  | Phi weights       | 386.4   | 524.9    | 530        | 532.1   | 533.3  |
  --------------------------------------------------------------------------

  2604 reflections have been flagged as outliers

  Summary statistics for 39747 observations matched to predictions:
  -----------------------------------------------------------------------------
  |                   | Min      | Q1       | Med        | Q3       | Max     |
  -----------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1008  | -0.02292 | 0.007548   | 0.03216  | 0.1216  |
  | Yc - Yo (mm)      | -0.08829 | -0.01313 | 0.005374   | 0.02507  | 0.09138 |
  | Phic - Phio (deg) | -0.05139 | -0.01036 | -0.0004898 | 0.009936 | 0.05098 |
  | X weights         | 231.1    | 397.5    | 402.5      | 404.5    | 405.6   |
  | Y weights         | 231.2    | 396.7    | 402.3      | 404.5    | 405.6   |
  | Phi weights       | 399.5    | 525      | 530        | 532.1    | 533.3   |
  -----------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.041131 | 0.028769 | 0.014744 |
  | 1    | 8099 | 0.038326 | 0.028061 | 0.015152 |
  | 2    | 8099 | 0.038311 | 0.028057 | 0.015084 |
  | 3    | 8099 | 0.038287 | 0.028028 | 0.01499  |
  | 4    | 8099 | 0.038244 | 0.02794  | 0.014908 |
  | 5    | 8099 | 0.038169 | 0.027766 | 0.014871 |
  | 6    | 8099 | 0.038061 | 0.027586 | 0.014865 |
  | 7    | 8099 | 0.037957 | 0.027548 | 0.014872 |
  | 8    | 8099 | 0.037909 | 0.027583 | 0.014881 |
  | 9    | 8099 | 0.0379   | 0.027593 | 0.014883 |
  | 10   | 8099 | 0.0379   | 0.027593 | 0.014883 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 8099 | 0.22035 | 0.16043 | 0.099218 |
  ---------------------------------------------
  Increasing resolution to 1.9 Angstrom
  model 1 (89300 reflections):
  Crystal:
      Unit cell: (57.816, 57.781, 150.030, 90.020, 89.997, 89.987)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2588, -0.9020},
                  { 0.8914,  0.3909,  0.2293},
                  { 0.2932, -0.8833,  0.3658}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  {-0.0000,  0.0000,  0.0067}}
      A = UB:    {{ 0.0060, -0.0045, -0.0060},
                  { 0.0154,  0.0068,  0.0015},
                  { 0.0051, -0.0153,  0.0024}}


  300 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 4)
  ################################################################################


  Summary statistics for 88041 observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.4767 | -0.02958 | 0.001571   | 0.03312 | 0.6737 |
  | Yc - Yo (mm)      | -1.42   | -0.01996 | 0.003176   | 0.02597 | 1.432  |
  | Phic - Phio (deg) | -1.434  | -0.01299 | -0.0004573 | 0.01265 | 0.9043 |
  | X weights         | 202.8   | 384.6    | 397.4      | 403.3   | 405.6  |
  | Y weights         | 171     | 378.6    | 395        | 402.6   | 405.6  |
  | Phi weights       | 318.9   | 520.3    | 529.4      | 533.3   | 533.3  |
  --------------------------------------------------------------------------

  8384 reflections have been flagged as outliers

  Summary statistics for 79657 observations matched to predictions:
  ---------------------------------------------------------------------------
  |                   | Min      | Q1       | Med       | Q3      | Max     |
  ---------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1483  | -0.028   | 0.002817  | 0.03259 | 0.1453  |
  | Yc - Yo (mm)      | -0.09244 | -0.01713 | 0.003795  | 0.02505 | 0.09839 |
  | Phic - Phio (deg) | -0.05483 | -0.01184 | -0.000235 | 0.01191 | 0.0556  |
  | X weights         | 205      | 387.5    | 398.3     | 403.5   | 405.6   |
  | Y weights         | 171      | 382.7    | 396.4     | 402.9   | 405.6   |
  | Phi weights       | 318.9    | 520.4    | 529.1     | 533     | 533.3   |
  ---------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.042412 | 0.033194 | 0.017463 |
  | 1    | 8099 | 0.041445 | 0.032572 | 0.017543 |
  | 2    | 8099 | 0.041436 | 0.032546 | 0.017495 |
  | 3    | 8099 | 0.041429 | 0.032461 | 0.017445 |
  | 4    | 8099 | 0.041425 | 0.032271 | 0.017404 |
  | 5    | 8099 | 0.041458 | 0.031963 | 0.017375 |
  | 6    | 8099 | 0.041544 | 0.031644 | 0.017351 |
  | 7    | 8099 | 0.041621 | 0.031457 | 0.017343 |
  | 8    | 8099 | 0.041651 | 0.031405 | 0.017343 |
  | 9    | 8099 | 0.041656 | 0.031398 | 0.017343 |
  | 10   | 8099 | 0.041656 | 0.031398 | 0.017343 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 8099 | 0.24219 | 0.18254 | 0.11562  |
  ---------------------------------------------
  Increasing resolution to 1.3 Angstrom
  model 1 (114690 reflections):
  Crystal:
      Unit cell: (57.815, 57.783, 150.035, 90.016, 89.991, 89.988)
      Space group: P 1
      U matrix:  {{ 0.3455, -0.2589, -0.9020},
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
  Starting refinement (macro-cycle 5)
  ################################################################################


  Summary statistics for 113247 observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5332 | -0.03414 | -0.004714  | 0.02995 | 0.6695 |
  | Yc - Yo (mm)      | -1.422  | -0.03176 | -0.003558  | 0.01999 | 1.273  |
  | Phic - Phio (deg) | -1.427  | -0.01451 | -0.0005735 | 0.01382 | 0.9041 |
  | X weights         | 135.2   | 371.6    | 393.4      | 402.6   | 405.6  |
  | Y weights         | 153.1   | 361.8    | 389        | 401.2   | 405.6  |
  | Phi weights       | 318.9   | 519.5    | 530.4      | 533.3   | 533.3  |
  --------------------------------------------------------------------------

  15804 reflections have been flagged as outliers

  Summary statistics for 97443 observations matched to predictions:
  ----------------------------------------------------------------------------
  |                   | Min      | Q1       | Med        | Q3      | Max     |
  ----------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1543  | -0.03126 | -0.003093  | 0.02881 | 0.1455  |
  | Yc - Yo (mm)      | -0.1028  | -0.02321 | -0.0006546 | 0.0205  | 0.1025  |
  | Phic - Phio (deg) | -0.06064 | -0.01269 | -0.0003442 | 0.01256 | 0.06094 |
  | X weights         | 135.2    | 379.3    | 395.7      | 403     | 405.6   |
  | Y weights         | 162.2    | 372.1    | 392.6      | 402     | 405.6   |
  | Phi weights       | 318.9    | 519.4    | 529.7      | 533.3   | 533.3   |
  ----------------------------------------------------------------------------


  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.043616 | 0.033904 | 0.019225 |
  | 1    | 8099 | 0.043394 | 0.033843 | 0.019184 |
  | 2    | 8099 | 0.043371 | 0.033778 | 0.019175 |
  | 3    | 8099 | 0.043321 | 0.033677 | 0.019169 |
  | 4    | 8099 | 0.043235 | 0.033557 | 0.019171 |
  | 5    | 8099 | 0.043138 | 0.033443 | 0.019179 |
  | 6    | 8099 | 0.043081 | 0.033365 | 0.019197 |
  | 7    | 8099 | 0.043068 | 0.033318 | 0.019215 |
  | 8    | 8099 | 0.043069 | 0.033302 | 0.019224 |
  | 9    | 8099 | 0.04307  | 0.0333   | 0.019225 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  --------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y | RMSD_Z   |
  |     |      | (px)    | (px)   | (images) |
  --------------------------------------------
  | 0   | 8099 | 0.25041 | 0.1936 | 0.12817  |
  --------------------------------------------
  Final refined crystal models:
  model 1 (114690 reflections):
  Crystal:
      Unit cell: (57.814, 57.784, 150.035, 90.012, 89.989, 89.988)
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

  Saving refined experiments to experiments.json
  Saving refined reflections to indexed.pickle


It is worth looking through this output to understand what the indexing program
has done. Note that this log
is automatically captured in the file :file:`dials.index.log`. There is also
a somewhat more information written into :file:`dials.index.debug.log`, but
this is probably only helpful if something has gone wrong and you are trying
to track down why.

Inspecting the log shows that the indexing step is done at fairly low
resolution: ``Setting d_min: 3.89``. The resolution limit of data that
can be used in indexing is determined by the size of the 3D FFT grid and the
likely maximum cell dimension. Here we
used :math:`256^3` grid points: ``FFT gridding: (256,256,256)``.
What follows are four macrocycles
of refinement at increasing resolution to bootstrap the indexing solution to as
many of the strong reflections as possible. In each case you can see that only
8099 reflections are used in the refinement job. The diffraction geometry is
here described by only 16 parameters (6 for the detector, 1 beam angle, 3
crystal 'misset' angles and 6 triclinic cell parameters). The problem is thus
hugely overdetermined. In order to save time, refinement uses a subset of the
input reflections, by default using 100 reflections for every degree of the scan.

Continuing to look through the log, we see that the first macrocyle of refinement makes
a big improvement in the positional RMSDs. The second macrocycle includes more reflections, after
extending to 3.2 Angstroms. The current model now shows slightly worse RMSDs
at the start, now that the higher resolution reflections are included, but refinement reduces
these again.
A similar situation is observed on the third and fourth macrocycles.
The RMSDs start higher again, now that more reflections are included, but refinement
is able to drive these down a little.
The final macrocycle includes data out to 1.3 Angstroms and refinement produces
a final model with
RMSDs of 0.043 mm in X, 0.033 mm in Y and 0.019 degrees in :math:`\phi`, corresponding
to 0.25 pixels in X, 0.19 pixels in Y and 0.13 image widths in :math:`\phi`.

Despite the high quality of this data, we notice from the ``Summary statistics``
tables that there were some outliers identified and removed from
refinement as resolution increases.
In the final macrocyle, prior to outlier rejection, we see the
distribution of positional residuals in the Y direction is tight around the
median, except for extreme values both positive and negative of more than 1 mm.
The angular residuals show a similar pattern with half the data having residuals
of less than about 0.14 degrees from the predicted positions, but the extreme
is as much as 1.4 degrees from the predicted diffraction angle. Large outliers
can dominate refinement using a least squares target, so it is important
to be able to remove these.

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

gives a table containing scoring data and unit cell for
each Bravais setting. The scores include the the metric fit (in degrees),
RMSDs (in mm), and the best and worse correlation coefficients for data
related by symmetry elements implied by the lowest symmetry space group from the
Bravais setting. This uses the raw spot intensity measurement from the
spot-finding procedure (uncorrected and unscaled) but provides a very
useful check to see if the data does appear to adhere to the proposed
symmetry operators.

::

  The following parameters have been modified:

  input {
    experiments = experiments.json
    reflections = indexed.pickle
  }

  -----------------------------------------------------------------------------------------------------------------
  Solution Metric fit  rmsd  min/max cc #spots lattice                                 unit_cell  volume      cb_op
  -----------------------------------------------------------------------------------------------------------------
         9     0.0336 0.060 0.787/0.848   8099      tP  57.79  57.79 150.01  90.00  90.00  90.00  500948      a,b,c
         8     0.0336 0.060 0.787/0.970   8099      oC  81.72  81.74 150.02  90.00  90.00  90.00 1002068  a-b,a+b,c
         7     0.0294 0.059 0.970/0.970   8099      mC  81.73  81.75 150.04  90.00  89.99  90.00 1002461  a-b,a+b,c
         6     0.0336 0.059 0.795/0.795   8099      mC  81.74  81.72 150.02  90.00  89.99  90.00 1002088 a+b,-a+b,c
         5     0.0167 0.059 0.787/0.899   8099      oP  57.80  57.77 150.01  90.00  90.00  90.00  500871      a,b,c
         4     0.0161 0.057 0.807/0.807   8099      mP  57.77  57.80 150.02  90.00  90.02  90.00  500971   -b,-a,-c
         3     0.0167 0.056 0.899/0.899   8099      mP  57.80  57.78 150.02  90.00  89.98  90.00  501005      a,b,c
         2     0.0163 0.057 0.787/0.787   8099      mP  57.78 150.01  57.80  90.00  89.99  90.00  501004      b,c,a
         1     0.0000 0.056         -/-   8099      aP  57.81  57.78 150.03  90.01  89.99  89.99  501195      a,b,c
  -----------------------------------------------------------------------------------------------------------------
  usr+sys time: 0.74 seconds
  wall clock time: 10.48 seconds


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
step using :doc:`dials.refine </programs/dials_refine>` in here. There
are many options to refinement. As an
aside, to show all the options up to and including ``expert_level = 1``
use this command::

  dials.refine -c -e 1

Equivalent command-line options exist for all the main DIALS programs.

The main reason
we may want to do an additional refinement job is to use a more sophisticated
model for the crystal,
allowing small misset rotations to occur over the course of the scan.
There are usually even small changes to the
cell dimensions (typically resulting in a net increase in cell volume) caused
by exposure to radiation during data collection. To account for both of these
effects we can extend our parameterisation to obtain a smoothed 'scan-varying'
model for both the crystal orientation and unit cell. To do this, we run a
further refinement job starting from the output of the previous job::

  dials.refine bravais_setting_9.json indexed.pickle scan_varying=true

The output for this job is::

  The following parameters have been modified:

  refinement {
    parameterisation {
      crystal {
        scan_varying = True
      }
    }
  }
  input {
    experiments = bravais_setting_9.json
    reflections = indexed.pickle
  }

  Configuring refiner

  Summary statistics for 113247 observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5257 | -0.03392 | -0.002959  | 0.03191 | 0.6438 |
  | Yc - Yo (mm)      | -1.402  | -0.03054 | -0.001345  | 0.02846 | 1.257  |
  | Phic - Phio (deg) | -1.409  | -0.01535 | -0.0008392 | 0.01382 | 0.9103 |
  | X weights         | 135.2   | 371.6    | 393.4      | 402.6   | 405.6  |
  | Y weights         | 153.1   | 361.8    | 389        | 401.2   | 405.6  |
  | Phi weights       | 318.9   | 519.5    | 530.4      | 533.3   | 533.3  |
  --------------------------------------------------------------------------

  11434 reflections have been flagged as outliers

  Summary statistics for 101813 observations matched to predictions:
  ----------------------------------------------------------------------------
  |                   | Min      | Q1       | Med        | Q3      | Max     |
  ----------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.1536  | -0.03252 | -0.002648  | 0.03021 | 0.1385  |
  | Yc - Yo (mm)      | -0.1259  | -0.02752 | -0.0008168 | 0.02679 | 0.1165  |
  | Phic - Phio (deg) | -0.05984 | -0.01327 | -0.0003616 | 0.01312 | 0.06268 |
  | X weights         | 135.2    | 377.6    | 395.1      | 402.9   | 405.6   |
  | Y weights         | 153.1    | 369.3    | 391.6      | 401.8   | 405.6   |
  | Phi weights       | 318.9    | 519.2    | 529.8      | 533.3   | 533.3   |
  ----------------------------------------------------------------------------

  Performing refinement...

  Refinement steps:
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 8099 | 0.045699 | 0.039879 | 0.019944 |
  | 1    | 8099 | 0.04521  | 0.038045 | 0.019945 |
  | 2    | 8099 | 0.045239 | 0.037898 | 0.019908 |
  | 3    | 8099 | 0.045257 | 0.037767 | 0.01988  |
  | 4    | 8099 | 0.045264 | 0.037683 | 0.019781 |
  | 5    | 8099 | 0.045269 | 0.037656 | 0.019627 |
  | 6    | 8099 | 0.045271 | 0.037653 | 0.019546 |
  | 7    | 8099 | 0.045271 | 0.037653 | 0.019535 |
  | 8    | 8099 | 0.045271 | 0.037652 | 0.019534 |
  ------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment:
  --------------------------------------------
  | Exp | Nref | RMSD_X | RMSD_Y  | RMSD_Z   |
  |     |      | (px)   | (px)    | (images) |
  --------------------------------------------
  | 0   | 8099 | 0.2632 | 0.21891 | 0.13022  |
  --------------------------------------------
  Saving refined experiments to refined_experiments.json
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

Diffraction geometry refinement, even with a scan-varying crystal model,
is hugely over-determined. So it is reasonable
to use a small subset of the total number of reflections to refine the
model. However, if we are being extra careful about data processing
and don't mind a slightly longer run time we might want to use all reflections
instead. In that case, we could use the following command::

  dials.refine refined_experiments.json indexed.pickle scan_varying=true use_all_reflections=true

This improves on the positional RMSDs from the previous job, but the angular
RMSD is slightly worse. In any case, the differences are only in the third
decimal place::

  RMSDs by experiment:
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 101878 | 0.26205 | 0.21742 | 0.13255  |
  -----------------------------------------------

The actual effect on the integrated data (in this case) of using
the model refined against the full set of strong spots rather than
the 100 reflections per degree subset quality is negligible.

To view the smoothly varying crystal cell parameters use the following command::

  dials.plot_scan_varying_crystal refined_experiments.json

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

  dials.integrate refined_experiments.json refined.pickle \
  outlier.algorithm=null nproc=4

The log file is quite long.

.. container:: toggle

    .. container:: header

        **Show/Hide Log**

    ::

      The following parameters have been modified:

      integration {
        mp {
          nproc = 4
        }
        background {
          simple {
            outlier {
              algorithm = *null nsigma truncated normal mosflm tukey
            }
          }
        }
      }
      input {
        experiments = refined_experiments.json
        reflections = refined.pickle
      }

      ================================================================================

      Initialising

      Processing reference reflections
       read 114690 strong spots
       using 114690 indexed reflections
       time taken: 0.0376


      ================================================================================

      Predicting reflections

      Prediction type: scan varying prediction
      Predicted 374050 reflections
      Matching reference spots with predicted reflections
       114690 observed reflections input
       374050 reflections predicted
       114526 reflections matched
       114525 reflections accepted
      Calculating E.S.D Beam Divergence.
      Calculating E.S.D Reflecting Range.
       sigma b: 0.021875 degrees
       sigma m: 0.065886 degrees

      ================================================================================

      Processing reflections

       Processing the following experiments:

       Experiments: 1
       Beams:       1
       Detectors:   1
       Goniometers: 1
       Scans:       1
       Crystals:    1
       Imagesets:   1

      ================================================================================

      Modelling reflection profiles

       Split 1254 reflections overlapping job boundaries

      Processing reflections in the following blocks of images:

       block_size: 42 frames

       --------------------------------------------------------------------------
        # | Group | Frame From | Frame To | Angle From | Angle To | # Reflections
       --------------------------------------------------------------------------
        0 |     0 |          0 |       42 |       82.0 |     88.3 |          6737
        1 |     0 |         21 |       63 |      85.15 |    91.45 |          4389
        2 |     0 |         42 |       84 |       88.3 |     94.6 |          4455
        3 |     0 |         63 |      105 |      91.45 |    97.75 |          4375
        4 |     0 |         84 |      126 |       94.6 |    100.9 |          4463
        5 |     0 |        105 |      147 |      97.75 |   104.05 |          4430
        6 |     0 |        126 |      168 |      100.9 |    107.2 |          4396
        7 |     0 |        147 |      189 |     104.05 |   110.35 |          4503
        8 |     0 |        168 |      210 |      107.2 |    113.5 |          4430
        9 |     0 |        189 |      231 |     110.35 |   116.65 |          4435
       10 |     0 |        210 |      252 |      113.5 |    119.8 |          4428
       11 |     0 |        231 |      273 |     116.65 |   122.95 |          4443
       12 |     0 |        252 |      294 |      119.8 |    126.1 |          4461
       13 |     0 |        273 |      315 |     122.95 |   129.25 |          4549
       14 |     0 |        294 |      336 |      126.1 |    132.4 |          4499
       15 |     0 |        315 |      357 |     129.25 |   135.55 |          4504
       16 |     0 |        336 |      378 |      132.4 |    138.7 |          4568
       17 |     0 |        357 |      399 |     135.55 |   141.85 |          4557
       18 |     0 |        378 |      420 |      138.7 |    145.0 |          4618
       19 |     0 |        399 |      441 |     141.85 |   148.15 |          4529
       20 |     0 |        420 |      462 |      145.0 |    151.3 |          4471
       21 |     0 |        441 |      483 |     148.15 |   154.45 |          4511
       22 |     0 |        462 |      504 |      151.3 |    157.6 |          4443
       23 |     0 |        483 |      525 |     154.45 |   160.75 |          4196
       24 |     0 |        504 |      540 |      157.6 |    163.0 |          5680
       --------------------------------------------------------------------------

       Using multiprocessing with 4 parallel job(s) and 1 thread(s) per job

       Beginning modelling job 0

       Frames: 0 -> 42

       Number of reflections
        Partial:     432
        Full:        6305
        In ice ring: 0
        Total:       6737

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       0  [635 ]: ***********************************
       1  [864 ]: ************************************************
       2  [1006]: ********************************************************
       3  [1028]: *********************************************************
       4  [1012]: ********************************************************
       5  [1018]: *********************************************************
       6  [1009]: ********************************************************
       7  [1024]: *********************************************************
       8  [1071]: ************************************************************
       9  [1057]: ***********************************************************
       10 [1085]: ************************************************************
       11 [1082]: ************************************************************
       12 [1063]: ***********************************************************
       13 [1072]: ************************************************************
       14 [1091]: *************************************************************
       15 [1112]: **************************************************************
       16 [1143]: ****************************************************************
       17 [1160]: *****************************************************************
       18 [1166]: *****************************************************************
       19 [1213]: ********************************************************************
       20 [1183]: ******************************************************************
       21 [1169]: *****************************************************************
       22 [1150]: ****************************************************************
       23 [1176]: *****************************************************************
       24 [1154]: ****************************************************************
       25 [1175]: *****************************************************************
       26 [1146]: ****************************************************************
       27 [1144]: ****************************************************************
       28 [1119]: **************************************************************
       29 [1091]: *************************************************************
       30 [879 ]: *************************************************
       31 [642 ]: ***********************************
       32 [423 ]: ***********************
       33 [205 ]: ***********
       34 [117 ]: ******
       35 [74  ]: ****
       36 [53  ]: **
       37 [41  ]: **
       38 [33  ]: *
       39 [25  ]: *
       40 [22  ]: *
       41 [21  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0136942 GB

       Modelled     0 /    65 reflection profiles on image 1
       Modelled    52 /   177 reflection profiles on image 2
       Modelled   128 /   182 reflection profiles on image 3
       Modelled   177 /   203 reflection profiles on image 4
       Modelled   179 /   196 reflection profiles on image 5
       Modelled   177 /   191 reflection profiles on image 6
       Modelled   158 /   175 reflection profiles on image 7
       Modelled   186 /   206 reflection profiles on image 8
       Modelled   171 /   178 reflection profiles on image 9
       Modelled   206 /   218 reflection profiles on image 10
       Modelled   217 /   230 reflection profiles on image 11
       Modelled   200 /   210 reflection profiles on image 12
       Modelled   187 /   200 reflection profiles on image 13
       Modelled   191 /   202 reflection profiles on image 14
       Modelled   186 /   195 reflection profiles on image 15
       Modelled   202 /   210 reflection profiles on image 16
       Modelled   194 /   201 reflection profiles on image 17
       Modelled   172 /   181 reflection profiles on image 18
       Modelled   244 /   254 reflection profiles on image 19
       Modelled   224 /   231 reflection profiles on image 20
       Modelled   215 /   226 reflection profiles on image 21
       Modelled   184 /   189 reflection profiles on image 22
       Modelled   223 /   233 reflection profiles on image 23
       Modelled   191 /   198 reflection profiles on image 24
       Modelled   239 /   248 reflection profiles on image 25
       Modelled   204 /   215 reflection profiles on image 26
       Modelled   217 /   221 reflection profiles on image 27
       Modelled   182 /   192 reflection profiles on image 28
       Modelled   222 /   231 reflection profiles on image 29
       Modelled   227 /   237 reflection profiles on image 30
       Modelled   208 /   219 reflection profiles on image 31
       Modelled   211 /   218 reflection profiles on image 32
       Modelled    83 /    88 reflection profiles on image 33
       Modelled    42 /    43 reflection profiles on image 34
       Modelled    20 /    21 reflection profiles on image 35
       Modelled    12 /    12 reflection profiles on image 36
       Modelled     7 /     8 reflection profiles on image 37
       Modelled     6 /     8 reflection profiles on image 38
       Modelled     3 /     3 reflection profiles on image 39
       Modelled     1 /     1 reflection profiles on image 40
       Modelled    11 /    21 reflection profiles on image 41
       Beginning modelling job 1

       Frames: 21 -> 63

       Number of reflections
        Partial:     45
        Full:        4344
        In ice ring: 0
        Total:       4389

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       21 [1   ]:
       22 [8   ]:
       23 [10  ]:
       24 [12  ]:
       25 [17  ]: *
       26 [24  ]: *
       27 [36  ]: **
       28 [59  ]: ***
       29 [94  ]: *****
       30 [267 ]: ***************
       31 [460 ]: ***************************
       32 [654 ]: **************************************
       33 [858 ]: **************************************************
       34 [947 ]: *******************************************************
       35 [986 ]: **********************************************************
       36 [1025]: ************************************************************
       37 [1044]: *************************************************************
       38 [1054]: **************************************************************
       39 [1034]: ************************************************************
       40 [1021]: ************************************************************
       41 [1053]: *************************************************************
       42 [1062]: **************************************************************
       43 [1118]: *****************************************************************
       44 [1124]: ******************************************************************
       45 [1156]: ********************************************************************
       46 [1127]: ******************************************************************
       47 [1110]: *****************************************************************
       48 [1096]: ****************************************************************
       49 [1105]: *****************************************************************
       50 [1027]: ************************************************************
       51 [822 ]: ************************************************
       52 [599 ]: ***********************************
       53 [355 ]: ********************
       54 [173 ]: **********
       55 [92  ]: *****
       56 [59  ]: ***
       57 [45  ]: **
       58 [38  ]: **
       59 [30  ]: *
       60 [25  ]: *
       61 [22  ]: *
       62 [19  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0122074 GB

       Modelled   105 /   108 reflection profiles on image 33
       Modelled   154 /   159 reflection profiles on image 34
       Modelled   170 /   182 reflection profiles on image 35
       Modelled   168 /   174 reflection profiles on image 36
       Modelled   187 /   194 reflection profiles on image 37
       Modelled   194 /   199 reflection profiles on image 38
       Modelled   212 /   223 reflection profiles on image 39
       Modelled   184 /   190 reflection profiles on image 40
       Modelled   220 /   240 reflection profiles on image 41
       Modelled   169 /   178 reflection profiles on image 42
       Modelled   198 /   207 reflection profiles on image 43
       Modelled   174 /   178 reflection profiles on image 44
       Modelled   215 /   221 reflection profiles on image 45
       Modelled   230 /   240 reflection profiles on image 46
       Modelled   213 /   219 reflection profiles on image 47
       Modelled   204 /   211 reflection profiles on image 48
       Modelled   217 /   225 reflection profiles on image 49
       Modelled   212 /   219 reflection profiles on image 50
       Modelled   212 /   223 reflection profiles on image 51
       Modelled   235 /   244 reflection profiles on image 52
       Modelled   174 /   182 reflection profiles on image 53
       Modelled    76 /    81 reflection profiles on image 54
       Modelled    29 /    33 reflection profiles on image 55
       Modelled    12 /    14 reflection profiles on image 56
       Modelled     7 /     7 reflection profiles on image 57
       Modelled     8 /     8 reflection profiles on image 58
       Modelled     4 /     5 reflection profiles on image 59
       Modelled     3 /     3 reflection profiles on image 60
       Modelled     3 /     3 reflection profiles on image 61
       Modelled     2 /    19 reflection profiles on image 62
       Beginning modelling job 2

       Frames: 42 -> 84

       Number of reflections
        Partial:     32
        Full:        4423
        In ice ring: 0
        Total:       4455

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       42 [3   ]:
       43 [7   ]:
       44 [9   ]:
       45 [18  ]: *
       46 [24  ]: *
       47 [33  ]: *
       48 [44  ]: **
       49 [61  ]: ***
       50 [103 ]: ******
       51 [289 ]: *****************
       52 [493 ]: *****************************
       53 [688 ]: *****************************************
       54 [901 ]: ******************************************************
       55 [981 ]: ***********************************************************
       56 [1048]: ***************************************************************
       57 [1066]: ****************************************************************
       58 [1101]: ******************************************************************
       59 [1113]: *******************************************************************
       60 [1103]: ******************************************************************
       61 [1090]: *****************************************************************
       62 [1047]: ***************************************************************
       63 [1053]: ***************************************************************
       64 [1092]: *****************************************************************
       65 [1110]: ******************************************************************
       66 [1127]: ********************************************************************
       67 [1114]: *******************************************************************
       68 [1106]: ******************************************************************
       69 [1064]: ****************************************************************
       70 [1049]: ***************************************************************
       71 [1019]: *************************************************************
       72 [834 ]: **************************************************
       73 [603 ]: ************************************
       74 [402 ]: ************************
       75 [163 ]: *********
       76 [83  ]: *****
       77 [59  ]: ***
       78 [39  ]: **
       79 [30  ]: *
       80 [22  ]: *
       81 [17  ]: *
       82 [11  ]:
       83 [9   ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0119831 GB

       Modelled   118 /   125 reflection profiles on image 54
       Modelled   151 /   159 reflection profiles on image 55
       Modelled   194 /   199 reflection profiles on image 56
       Modelled   190 /   198 reflection profiles on image 57
       Modelled   190 /   194 reflection profiles on image 58
       Modelled   202 /   205 reflection profiles on image 59
       Modelled   208 /   218 reflection profiles on image 60
       Modelled   231 /   242 reflection profiles on image 61
       Modelled   205 /   220 reflection profiles on image 62
       Modelled   176 /   186 reflection profiles on image 63
       Modelled   179 /   184 reflection profiles on image 64
       Modelled   192 /   197 reflection profiles on image 65
       Modelled   197 /   207 reflection profiles on image 66
       Modelled   214 /   223 reflection profiles on image 67
       Modelled   221 /   231 reflection profiles on image 68
       Modelled   202 /   206 reflection profiles on image 69
       Modelled   193 /   206 reflection profiles on image 70
       Modelled   214 /   221 reflection profiles on image 71
       Modelled   224 /   231 reflection profiles on image 72
       Modelled   194 /   201 reflection profiles on image 73
       Modelled   229 /   239 reflection profiles on image 74
       Modelled    77 /    80 reflection profiles on image 75
       Modelled    23 /    24 reflection profiles on image 76
       Modelled    19 /    20 reflection profiles on image 77
       Modelled     7 /     9 reflection profiles on image 78
       Modelled     7 /     8 reflection profiles on image 79
       Modelled     5 /     5 reflection profiles on image 80
       Modelled     4 /     6 reflection profiles on image 81
       Modelled     2 /     2 reflection profiles on image 82
       Modelled     0 /     9 reflection profiles on image 83
       Beginning modelling job 3

       Frames: 63 -> 105

       Number of reflections
        Partial:     44
        Full:        4331
        In ice ring: 0
        Total:       4375

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       63  [3   ]:
       64  [7   ]:
       65  [11  ]:
       66  [15  ]:
       67  [24  ]: *
       68  [29  ]: *
       69  [47  ]: **
       70  [79  ]: ****
       71  [139 ]: *******
       72  [334 ]: ******************
       73  [526 ]: *****************************
       74  [707 ]: ****************************************
       75  [900 ]: ***************************************************
       76  [955 ]: ******************************************************
       77  [979 ]: *******************************************************
       78  [1003]: ********************************************************
       79  [1022]: **********************************************************
       80  [1058]: ************************************************************
       81  [1088]: *************************************************************
       82  [1102]: **************************************************************
       83  [1140]: ****************************************************************
       84  [1167]: ******************************************************************
       85  [1180]: *******************************************************************
       86  [1178]: ******************************************************************
       87  [1143]: ****************************************************************
       88  [1074]: ************************************************************
       89  [1038]: **********************************************************
       90  [1016]: *********************************************************
       91  [1010]: *********************************************************
       92  [948 ]: *****************************************************
       93  [808 ]: *********************************************
       94  [605 ]: **********************************
       95  [388 ]: **********************
       96  [188 ]: **********
       97  [103 ]: *****
       98  [66  ]: ***
       99  [51  ]: **
       100 [39  ]: **
       101 [34  ]: *
       102 [27  ]: *
       103 [22  ]: *
       104 [19  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0127539 GB

       Modelled   120 /   124 reflection profiles on image 75
       Modelled   162 /   168 reflection profiles on image 76
       Modelled   177 /   182 reflection profiles on image 77
       Modelled   185 /   192 reflection profiles on image 78
       Modelled   178 /   183 reflection profiles on image 79
       Modelled   188 /   195 reflection profiles on image 80
       Modelled   180 /   190 reflection profiles on image 81
       Modelled   192 /   206 reflection profiles on image 82
       Modelled   192 /   217 reflection profiles on image 83
       Modelled   206 /   215 reflection profiles on image 84
       Modelled   207 /   212 reflection profiles on image 85
       Modelled   212 /   220 reflection profiles on image 86
       Modelled   246 /   253 reflection profiles on image 87
       Modelled   202 /   211 reflection profiles on image 88
       Modelled   217 /   222 reflection profiles on image 89
       Modelled   180 /   192 reflection profiles on image 90
       Modelled   211 /   219 reflection profiles on image 91
       Modelled   159 /   166 reflection profiles on image 92
       Modelled   192 /   203 reflection profiles on image 93
       Modelled   215 /   217 reflection profiles on image 94
       Modelled   189 /   200 reflection profiles on image 95
       Modelled    84 /    85 reflection profiles on image 96
       Modelled    37 /    37 reflection profiles on image 97
       Modelled    15 /    15 reflection profiles on image 98
       Modelled    12 /    12 reflection profiles on image 99
       Modelled     5 /     5 reflection profiles on image 100
       Modelled     6 /     7 reflection profiles on image 101
       Modelled     5 /     5 reflection profiles on image 102
       Modelled     2 /     3 reflection profiles on image 103
       Modelled     1 /    19 reflection profiles on image 104
       Beginning modelling job 4

       Frames: 84 -> 126

       Number of reflections
        Partial:     46
        Full:        4417
        In ice ring: 0
        Total:       4463

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       84  [1   ]:
       85  [2   ]:
       86  [5   ]:
       87  [9   ]:
       88  [14  ]:
       89  [25  ]: *
       90  [36  ]: **
       91  [62  ]: ***
       92  [105 ]: ******
       93  [309 ]: ******************
       94  [530 ]: *******************************
       95  [741 ]: ********************************************
       96  [967 ]: *********************************************************
       97  [1070]: ****************************************************************
       98  [1086]: *****************************************************************
       99  [1092]: *****************************************************************
       100 [1110]: ******************************************************************
       101 [1071]: ****************************************************************
       102 [1076]: ****************************************************************
       103 [1068]: ***************************************************************
       104 [1082]: ****************************************************************
       105 [1119]: *******************************************************************
       106 [1118]: ******************************************************************
       107 [1079]: ****************************************************************
       108 [1063]: ***************************************************************
       109 [1060]: ***************************************************************
       110 [1069]: ****************************************************************
       111 [1059]: ***************************************************************
       112 [1042]: **************************************************************
       113 [1019]: *************************************************************
       114 [810 ]: ************************************************
       115 [602 ]: ************************************
       116 [408 ]: ************************
       117 [177 ]: **********
       118 [96  ]: *****
       119 [68  ]: ****
       120 [45  ]: **
       121 [33  ]: *
       122 [29  ]: *
       123 [27  ]: *
       124 [18  ]: *
       125 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0122721 GB

       Modelled   138 /   145 reflection profiles on image 96
       Modelled   174 /   179 reflection profiles on image 97
       Modelled   204 /   209 reflection profiles on image 98
       Modelled   178 /   180 reflection profiles on image 99
       Modelled   233 /   242 reflection profiles on image 100
       Modelled   197 /   207 reflection profiles on image 101
       Modelled   212 /   217 reflection profiles on image 102
       Modelled   185 /   191 reflection profiles on image 103
       Modelled   189 /   205 reflection profiles on image 104
       Modelled   186 /   195 reflection profiles on image 105
       Modelled   209 /   224 reflection profiles on image 106
       Modelled   207 /   214 reflection profiles on image 107
       Modelled   204 /   217 reflection profiles on image 108
       Modelled   192 /   200 reflection profiles on image 109
       Modelled   183 /   189 reflection profiles on image 110
       Modelled   190 /   202 reflection profiles on image 111
       Modelled   187 /   199 reflection profiles on image 112
       Modelled   231 /   238 reflection profiles on image 113
       Modelled   201 /   208 reflection profiles on image 114
       Modelled   191 /   194 reflection profiles on image 115
       Modelled   225 /   231 reflection profiles on image 116
       Modelled    77 /    81 reflection profiles on image 117
       Modelled    28 /    28 reflection profiles on image 118
       Modelled    22 /    23 reflection profiles on image 119
       Modelled    11 /    12 reflection profiles on image 120
       Modelled     4 /     4 reflection profiles on image 121
       Modelled     1 /     2 reflection profiles on image 122
       Modelled     9 /     9 reflection profiles on image 123
       Modelled     1 /     2 reflection profiles on image 124
       Modelled     2 /    16 reflection profiles on image 125
       Beginning modelling job 5

       Frames: 105 -> 147

       Number of reflections
        Partial:     59
        Full:        4371
        In ice ring: 0
        Total:       4430

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       105 [3   ]:
       106 [5   ]:
       107 [9   ]:
       108 [10  ]:
       109 [15  ]:
       110 [27  ]: *
       111 [42  ]: **
       112 [65  ]: ***
       113 [118 ]: *******
       114 [287 ]: *****************
       115 [515 ]: ******************************
       116 [749 ]: ********************************************
       117 [953 ]: *********************************************************
       118 [1030]: *************************************************************
       119 [1049]: **************************************************************
       120 [1037]: **************************************************************
       121 [1049]: **************************************************************
       122 [1068]: ****************************************************************
       123 [1084]: ****************************************************************
       124 [1056]: ***************************************************************
       125 [1058]: ***************************************************************
       126 [1112]: ******************************************************************
       127 [1076]: ****************************************************************
       128 [1118]: *******************************************************************
       129 [1113]: ******************************************************************
       130 [1110]: ******************************************************************
       131 [1097]: *****************************************************************
       132 [1082]: ****************************************************************
       133 [1035]: **************************************************************
       134 [968 ]: **********************************************************
       135 [781 ]: **********************************************
       136 [596 ]: ***********************************
       137 [365 ]: *********************
       138 [178 ]: **********
       139 [93  ]: *****
       140 [61  ]: ***
       141 [47  ]: **
       142 [44  ]: **
       143 [34  ]: **
       144 [28  ]: *
       145 [25  ]: *
       146 [24  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0119768 GB

       Modelled   112 /   114 reflection profiles on image 117
       Modelled   173 /   182 reflection profiles on image 118
       Modelled   199 /   204 reflection profiles on image 119
       Modelled   201 /   206 reflection profiles on image 120
       Modelled   192 /   195 reflection profiles on image 121
       Modelled   183 /   191 reflection profiles on image 122
       Modelled   206 /   212 reflection profiles on image 123
       Modelled   209 /   214 reflection profiles on image 124
       Modelled   172 /   203 reflection profiles on image 125
       Modelled   218 /   229 reflection profiles on image 126
       Modelled   169 /   176 reflection profiles on image 127
       Modelled   195 /   203 reflection profiles on image 128
       Modelled   226 /   232 reflection profiles on image 129
       Modelled   206 /   210 reflection profiles on image 130
       Modelled   205 /   217 reflection profiles on image 131
       Modelled   210 /   220 reflection profiles on image 132
       Modelled   221 /   227 reflection profiles on image 133
       Modelled   210 /   214 reflection profiles on image 134
       Modelled   175 /   185 reflection profiles on image 135
       Modelled   228 /   231 reflection profiles on image 136
       Modelled   179 /   187 reflection profiles on image 137
       Modelled    81 /    85 reflection profiles on image 138
       Modelled    31 /    32 reflection profiles on image 139
       Modelled    14 /    14 reflection profiles on image 140
       Modelled     3 /     3 reflection profiles on image 141
       Modelled     8 /    10 reflection profiles on image 142
       Modelled     6 /     6 reflection profiles on image 143
       Modelled     3 /     3 reflection profiles on image 144
       Modelled     1 /     1 reflection profiles on image 145
       Modelled     1 /    24 reflection profiles on image 146
       Beginning modelling job 6

       Frames: 126 -> 168

       Number of reflections
        Partial:     43
        Full:        4353
        In ice ring: 0
        Total:       4396

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       126 [1   ]:
       127 [2   ]:
       128 [3   ]:
       129 [4   ]:
       130 [13  ]:
       131 [23  ]: *
       132 [39  ]: **
       133 [69  ]: ****
       134 [127 ]: *******
       135 [350 ]: ********************
       136 [551 ]: ********************************
       137 [739 ]: *******************************************
       138 [925 ]: *******************************************************
       139 [1000]: ***********************************************************
       140 [1007]: ***********************************************************
       141 [1043]: **************************************************************
       142 [1060]: ***************************************************************
       143 [1044]: **************************************************************
       144 [1032]: *************************************************************
       145 [1061]: ***************************************************************
       146 [1069]: ***************************************************************
       147 [1106]: *****************************************************************
       148 [1122]: ******************************************************************
       149 [1112]: ******************************************************************
       150 [1119]: ******************************************************************
       151 [1126]: *******************************************************************
       152 [1119]: ******************************************************************
       153 [1080]: ****************************************************************
       154 [1052]: **************************************************************
       155 [995 ]: ***********************************************************
       156 [810 ]: ************************************************
       157 [578 ]: **********************************
       158 [395 ]: ***********************
       159 [176 ]: **********
       160 [97  ]: *****
       161 [64  ]: ***
       162 [47  ]: **
       163 [35  ]: **
       164 [25  ]: *
       165 [19  ]: *
       166 [16  ]:
       167 [14  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0120732 GB

       Modelled   132 /   136 reflection profiles on image 138
       Modelled   189 /   194 reflection profiles on image 139
       Modelled   147 /   152 reflection profiles on image 140
       Modelled   185 /   193 reflection profiles on image 141
       Modelled   204 /   208 reflection profiles on image 142
       Modelled   195 /   198 reflection profiles on image 143
       Modelled   180 /   186 reflection profiles on image 144
       Modelled   210 /   221 reflection profiles on image 145
       Modelled   181 /   196 reflection profiles on image 146
       Modelled   188 /   204 reflection profiles on image 147
       Modelled   203 /   211 reflection profiles on image 148
       Modelled   200 /   212 reflection profiles on image 149
       Modelled   198 /   201 reflection profiles on image 150
       Modelled   204 /   213 reflection profiles on image 151
       Modelled   214 /   223 reflection profiles on image 152
       Modelled   203 /   213 reflection profiles on image 153
       Modelled   207 /   217 reflection profiles on image 154
       Modelled   202 /   208 reflection profiles on image 155
       Modelled   227 /   232 reflection profiles on image 156
       Modelled   179 /   183 reflection profiles on image 157
       Modelled   210 /   219 reflection profiles on image 158
       Modelled    76 /    79 reflection profiles on image 159
       Modelled    32 /    33 reflection profiles on image 160
       Modelled    16 /    17 reflection profiles on image 161
       Modelled    11 /    12 reflection profiles on image 162
       Modelled     6 /    10 reflection profiles on image 163
       Modelled     6 /     6 reflection profiles on image 164
       Modelled     3 /     3 reflection profiles on image 165
       Modelled     1 /     2 reflection profiles on image 166
       Modelled     0 /    14 reflection profiles on image 167
       Beginning modelling job 7

       Frames: 147 -> 189

       Number of reflections
        Partial:     42
        Full:        4461
        In ice ring: 0
        Total:       4503

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       147 [3   ]:
       148 [8   ]:
       149 [12  ]:
       150 [18  ]: *
       151 [23  ]: *
       152 [31  ]: *
       153 [41  ]: **
       154 [66  ]: ***
       155 [119 ]: ******
       156 [330 ]: *******************
       157 [537 ]: *******************************
       158 [740 ]: *******************************************
       159 [952 ]: *******************************************************
       160 [1051]: *************************************************************
       161 [1085]: ***************************************************************
       162 [1087]: ***************************************************************
       163 [1097]: ****************************************************************
       164 [1121]: *****************************************************************
       165 [1120]: *****************************************************************
       166 [1146]: *******************************************************************
       167 [1142]: ******************************************************************
       168 [1129]: ******************************************************************
       169 [1124]: *****************************************************************
       170 [1134]: ******************************************************************
       171 [1138]: ******************************************************************
       172 [1109]: ****************************************************************
       173 [1079]: ***************************************************************
       174 [1044]: *************************************************************
       175 [1032]: ************************************************************
       176 [989 ]: *********************************************************
       177 [807 ]: ***********************************************
       178 [606 ]: ***********************************
       179 [390 ]: **********************
       180 [175 ]: **********
       181 [96  ]: *****
       182 [66  ]: ***
       183 [47  ]: **
       184 [38  ]: **
       185 [30  ]: *
       186 [26  ]: *
       187 [20  ]: *
       188 [18  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0125515 GB

       Modelled   131 /   135 reflection profiles on image 159
       Modelled   184 /   188 reflection profiles on image 160
       Modelled   196 /   203 reflection profiles on image 161
       Modelled   190 /   198 reflection profiles on image 162
       Modelled   202 /   210 reflection profiles on image 163
       Modelled   212 /   217 reflection profiles on image 164
       Modelled   191 /   200 reflection profiles on image 165
       Modelled   194 /   202 reflection profiles on image 166
       Modelled   209 /   233 reflection profiles on image 167
       Modelled   218 /   221 reflection profiles on image 168
       Modelled   196 /   206 reflection profiles on image 169
       Modelled   207 /   213 reflection profiles on image 170
       Modelled   210 /   219 reflection profiles on image 171
       Modelled   204 /   218 reflection profiles on image 172
       Modelled   216 /   222 reflection profiles on image 173
       Modelled   190 /   195 reflection profiles on image 174
       Modelled   207 /   212 reflection profiles on image 175
       Modelled   198 /   204 reflection profiles on image 176
       Modelled   189 /   201 reflection profiles on image 177
       Modelled   209 /   216 reflection profiles on image 178
       Modelled   211 /   215 reflection profiles on image 179
       Modelled    76 /    79 reflection profiles on image 180
       Modelled    30 /    30 reflection profiles on image 181
       Modelled    19 /    19 reflection profiles on image 182
       Modelled     9 /     9 reflection profiles on image 183
       Modelled     6 /     8 reflection profiles on image 184
       Modelled     3 /     4 reflection profiles on image 185
       Modelled     5 /     6 reflection profiles on image 186
       Modelled     2 /     2 reflection profiles on image 187
       Modelled     3 /    18 reflection profiles on image 188
       Beginning modelling job 8

       Frames: 168 -> 210

       Number of reflections
        Partial:     33
        Full:        4397
        In ice ring: 0
        Total:       4430

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       168 [1   ]:
       169 [6   ]:
       170 [9   ]:
       171 [10  ]:
       172 [13  ]:
       173 [18  ]: *
       174 [28  ]: *
       175 [55  ]: ***
       176 [108 ]: ******
       177 [288 ]: *****************
       178 [494 ]: *****************************
       179 [729 ]: ********************************************
       180 [943 ]: ********************************************************
       181 [1010]: ************************************************************
       182 [1038]: **************************************************************
       183 [1029]: **************************************************************
       184 [1046]: ***************************************************************
       185 [1073]: ****************************************************************
       186 [1067]: ****************************************************************
       187 [1083]: *****************************************************************
       188 [1094]: ******************************************************************
       189 [1101]: ******************************************************************
       190 [1095]: ******************************************************************
       191 [1081]: *****************************************************************
       192 [1089]: *****************************************************************
       193 [1110]: *******************************************************************
       194 [1104]: ******************************************************************
       195 [1091]: *****************************************************************
       196 [1060]: ***************************************************************
       197 [1008]: ************************************************************
       198 [798 ]: ************************************************
       199 [577 ]: **********************************
       200 [366 ]: **********************
       201 [178 ]: **********
       202 [82  ]: ****
       203 [56  ]: ***
       204 [37  ]: **
       205 [22  ]: *
       206 [18  ]: *
       207 [17  ]: *
       208 [14  ]:
       209 [13  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0116969 GB

       Modelled   129 /   133 reflection profiles on image 180
       Modelled   164 /   169 reflection profiles on image 181
       Modelled   183 /   190 reflection profiles on image 182
       Modelled   202 /   208 reflection profiles on image 183
       Modelled   185 /   189 reflection profiles on image 184
       Modelled   204 /   213 reflection profiles on image 185
       Modelled   187 /   189 reflection profiles on image 186
       Modelled   201 /   210 reflection profiles on image 187
       Modelled   199 /   217 reflection profiles on image 188
       Modelled   203 /   209 reflection profiles on image 189
       Modelled   218 /   228 reflection profiles on image 190
       Modelled   192 /   201 reflection profiles on image 191
       Modelled   197 /   203 reflection profiles on image 192
       Modelled   185 /   197 reflection profiles on image 193
       Modelled   220 /   223 reflection profiles on image 194
       Modelled   213 /   219 reflection profiles on image 195
       Modelled   196 /   200 reflection profiles on image 196
       Modelled   227 /   234 reflection profiles on image 197
       Modelled   213 /   221 reflection profiles on image 198
       Modelled   204 /   211 reflection profiles on image 199
       Modelled   179 /   188 reflection profiles on image 200
       Modelled    92 /    96 reflection profiles on image 201
       Modelled    26 /    26 reflection profiles on image 202
       Modelled    18 /    19 reflection profiles on image 203
       Modelled    15 /    15 reflection profiles on image 204
       Modelled     4 /     4 reflection profiles on image 205
       Modelled     1 /     1 reflection profiles on image 206
       Modelled     3 /     3 reflection profiles on image 207
       Modelled     1 /     1 reflection profiles on image 208
       Modelled     3 /    13 reflection profiles on image 209
       Beginning modelling job 9

       Frames: 189 -> 231

       Number of reflections
        Partial:     34
        Full:        4401
        In ice ring: 0
        Total:       4435

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       189 [3   ]:
       190 [4   ]:
       191 [6   ]:
       192 [8   ]:
       193 [12  ]:
       194 [21  ]: *
       195 [42  ]: **
       196 [61  ]: ***
       197 [117 ]: ******
       198 [282 ]: ****************
       199 [489 ]: ****************************
       200 [707 ]: *****************************************
       201 [920 ]: ******************************************************
       202 [1009]: ***********************************************************
       203 [1030]: ************************************************************
       204 [1048]: *************************************************************
       205 [1095]: ****************************************************************
       206 [1122]: ******************************************************************
       207 [1130]: ******************************************************************
       208 [1075]: ***************************************************************
       209 [1075]: ***************************************************************
       210 [1135]: *******************************************************************
       211 [1131]: ******************************************************************
       212 [1134]: ******************************************************************
       213 [1110]: *****************************************************************
       214 [1097]: ****************************************************************
       215 [1101]: ****************************************************************
       216 [1061]: **************************************************************
       217 [1022]: ************************************************************
       218 [994 ]: **********************************************************
       219 [790 ]: **********************************************
       220 [583 ]: **********************************
       221 [393 ]: ***********************
       222 [171 ]: **********
       223 [87  ]: *****
       224 [56  ]: ***
       225 [41  ]: **
       226 [33  ]: *
       227 [29  ]: *
       228 [26  ]: *
       229 [21  ]: *
       230 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0122668 GB

       Modelled   107 /   107 reflection profiles on image 201
       Modelled   180 /   185 reflection profiles on image 202
       Modelled   191 /   202 reflection profiles on image 203
       Modelled   184 /   192 reflection profiles on image 204
       Modelled   176 /   186 reflection profiles on image 205
       Modelled   197 /   202 reflection profiles on image 206
       Modelled   232 /   241 reflection profiles on image 207
       Modelled   197 /   205 reflection profiles on image 208
       Modelled   179 /   202 reflection profiles on image 209
       Modelled   203 /   213 reflection profiles on image 210
       Modelled   210 /   218 reflection profiles on image 211
       Modelled   201 /   208 reflection profiles on image 212
       Modelled   208 /   216 reflection profiles on image 213
       Modelled   216 /   219 reflection profiles on image 214
       Modelled   214 /   219 reflection profiles on image 215
       Modelled   205 /   214 reflection profiles on image 216
       Modelled   180 /   191 reflection profiles on image 217
       Modelled   220 /   225 reflection profiles on image 218
       Modelled   199 /   207 reflection profiles on image 219
       Modelled   178 /   190 reflection profiles on image 220
       Modelled   211 /   222 reflection profiles on image 221
       Modelled    82 /    84 reflection profiles on image 222
       Modelled    31 /    31 reflection profiles on image 223
       Modelled    13 /    15 reflection profiles on image 224
       Modelled     8 /     8 reflection profiles on image 225
       Modelled     3 /     4 reflection profiles on image 226
       Modelled     2 /     3 reflection profiles on image 227
       Modelled     4 /     5 reflection profiles on image 228
       Modelled     5 /     5 reflection profiles on image 229
       Modelled     3 /    16 reflection profiles on image 230
       Beginning modelling job 10

       Frames: 210 -> 252

       Number of reflections
        Partial:     45
        Full:        4383
        In ice ring: 0
        Total:       4428

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       210 [2   ]:
       211 [7   ]:
       212 [8   ]:
       213 [11  ]:
       214 [18  ]: *
       215 [25  ]: *
       216 [39  ]: **
       217 [60  ]: ***
       218 [103 ]: ******
       219 [281 ]: ****************
       220 [506 ]: ******************************
       221 [718 ]: *******************************************
       222 [944 ]: ********************************************************
       223 [1013]: *************************************************************
       224 [1066]: ****************************************************************
       225 [1081]: *****************************************************************
       226 [1103]: ******************************************************************
       227 [1092]: *****************************************************************
       228 [1110]: *******************************************************************
       229 [1097]: ******************************************************************
       230 [1058]: ***************************************************************
       231 [1081]: *****************************************************************
       232 [1051]: ***************************************************************
       233 [1028]: **************************************************************
       234 [1010]: ************************************************************
       235 [1031]: **************************************************************
       236 [1068]: ****************************************************************
       237 [1056]: ***************************************************************
       238 [1052]: ***************************************************************
       239 [1006]: ************************************************************
       240 [818 ]: *************************************************
       241 [613 ]: *************************************
       242 [418 ]: *************************
       243 [189 ]: ***********
       244 [102 ]: ******
       245 [66  ]: ***
       246 [50  ]: ***
       247 [36  ]: **
       248 [30  ]: *
       249 [27  ]: *
       250 [21  ]: *
       251 [20  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.011515 GB

       Modelled   126 /   129 reflection profiles on image 222
       Modelled   173 /   179 reflection profiles on image 223
       Modelled   190 /   196 reflection profiles on image 224
       Modelled   210 /   214 reflection profiles on image 225
       Modelled   188 /   195 reflection profiles on image 226
       Modelled   188 /   194 reflection profiles on image 227
       Modelled   226 /   235 reflection profiles on image 228
       Modelled   216 /   226 reflection profiles on image 229
       Modelled   188 /   215 reflection profiles on image 230
       Modelled   194 /   205 reflection profiles on image 231
       Modelled   205 /   215 reflection profiles on image 232
       Modelled   210 /   215 reflection profiles on image 233
       Modelled   187 /   194 reflection profiles on image 234
       Modelled   164 /   171 reflection profiles on image 235
       Modelled   183 /   192 reflection profiles on image 236
       Modelled   195 /   203 reflection profiles on image 237
       Modelled   208 /   220 reflection profiles on image 238
       Modelled   201 /   212 reflection profiles on image 239
       Modelled   196 /   205 reflection profiles on image 240
       Modelled   187 /   195 reflection profiles on image 241
       Modelled   223 /   229 reflection profiles on image 242
       Modelled    82 /    87 reflection profiles on image 243
       Modelled    35 /    36 reflection profiles on image 244
       Modelled    15 /    16 reflection profiles on image 245
       Modelled    13 /    14 reflection profiles on image 246
       Modelled     5 /     6 reflection profiles on image 247
       Modelled     3 /     3 reflection profiles on image 248
       Modelled     5 /     6 reflection profiles on image 249
       Modelled     0 /     1 reflection profiles on image 250
       Modelled     3 /    20 reflection profiles on image 251
       Beginning modelling job 11

       Frames: 231 -> 273

       Number of reflections
        Partial:     49
        Full:        4394
        In ice ring: 0
        Total:       4443

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       231 [1   ]:
       232 [1   ]:
       233 [4   ]:
       234 [15  ]:
       235 [20  ]: *
       236 [27  ]: *
       237 [44  ]: **
       238 [62  ]: ***
       239 [109 ]: ******
       240 [293 ]: *****************
       241 [510 ]: ******************************
       242 [725 ]: *******************************************
       243 [942 ]: ********************************************************
       244 [1020]: ************************************************************
       245 [1050]: **************************************************************
       246 [1062]: ***************************************************************
       247 [1092]: *****************************************************************
       248 [1106]: *****************************************************************
       249 [1097]: *****************************************************************
       250 [1097]: *****************************************************************
       251 [1099]: *****************************************************************
       252 [1102]: *****************************************************************
       253 [1125]: *******************************************************************
       254 [1111]: ******************************************************************
       255 [1059]: ***************************************************************
       256 [1037]: *************************************************************
       257 [1015]: ************************************************************
       258 [1046]: **************************************************************
       259 [1033]: *************************************************************
       260 [989 ]: **********************************************************
       261 [837 ]: *************************************************
       262 [615 ]: ************************************
       263 [411 ]: ************************
       264 [200 ]: ***********
       265 [101 ]: ******
       266 [67  ]: ***
       267 [49  ]: **
       268 [35  ]: **
       269 [30  ]: *
       270 [26  ]: *
       271 [25  ]: *
       272 [23  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0120038 GB

       Modelled   117 /   118 reflection profiles on image 243
       Modelled   177 /   184 reflection profiles on image 244
       Modelled   189 /   198 reflection profiles on image 245
       Modelled   181 /   189 reflection profiles on image 246
       Modelled   188 /   196 reflection profiles on image 247
       Modelled   217 /   223 reflection profiles on image 248
       Modelled   212 /   217 reflection profiles on image 249
       Modelled   216 /   224 reflection profiles on image 250
       Modelled   200 /   226 reflection profiles on image 251
       Modelled   183 /   192 reflection profiles on image 252
       Modelled   199 /   205 reflection profiles on image 253
       Modelled   232 /   238 reflection profiles on image 254
       Modelled   211 /   218 reflection profiles on image 255
       Modelled   195 /   203 reflection profiles on image 256
       Modelled   180 /   184 reflection profiles on image 257
       Modelled   200 /   205 reflection profiles on image 258
       Modelled   201 /   206 reflection profiles on image 259
       Modelled   170 /   180 reflection profiles on image 260
       Modelled   217 /   222 reflection profiles on image 261
       Modelled   198 /   204 reflection profiles on image 262
       Modelled   208 /   211 reflection profiles on image 263
       Modelled    93 /    99 reflection profiles on image 264
       Modelled    31 /    34 reflection profiles on image 265
       Modelled    17 /    18 reflection profiles on image 266
       Modelled    11 /    14 reflection profiles on image 267
       Modelled     5 /     5 reflection profiles on image 268
       Modelled     4 /     4 reflection profiles on image 269
       Modelled     1 /     1 reflection profiles on image 270
       Modelled     2 /     2 reflection profiles on image 271
       Modelled     5 /    23 reflection profiles on image 272
       Beginning modelling job 12

       Frames: 252 -> 294

       Number of reflections
        Partial:     35
        Full:        4426
        In ice ring: 0
        Total:       4461

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       252 [1   ]:
       253 [4   ]:
       254 [11  ]:
       255 [18  ]: *
       256 [23  ]: *
       257 [31  ]: *
       258 [46  ]: **
       259 [65  ]: ***
       260 [112 ]: ******
       261 [277 ]: ****************
       262 [489 ]: ****************************
       263 [698 ]: *****************************************
       264 [878 ]: ***************************************************
       265 [997 ]: **********************************************************
       266 [1032]: ************************************************************
       267 [1046]: *************************************************************
       268 [1089]: ****************************************************************
       269 [1057]: **************************************************************
       270 [1079]: ***************************************************************
       271 [1078]: ***************************************************************
       272 [1081]: ***************************************************************
       273 [1136]: *******************************************************************
       274 [1135]: ******************************************************************
       275 [1133]: ******************************************************************
       276 [1134]: ******************************************************************
       277 [1111]: *****************************************************************
       278 [1101]: ****************************************************************
       279 [1050]: *************************************************************
       280 [1018]: ************************************************************
       281 [1004]: ***********************************************************
       282 [836 ]: *************************************************
       283 [605 ]: ***********************************
       284 [405 ]: ***********************
       285 [193 ]: ***********
       286 [96  ]: *****
       287 [59  ]: ***
       288 [42  ]: **
       289 [27  ]: *
       290 [19  ]: *
       291 [15  ]:
       292 [15  ]:
       293 [14  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0120123 GB

       Modelled   108 /   113 reflection profiles on image 264
       Modelled   164 /   173 reflection profiles on image 265
       Modelled   186 /   189 reflection profiles on image 266
       Modelled   181 /   188 reflection profiles on image 267
       Modelled   205 /   209 reflection profiles on image 268
       Modelled   193 /   198 reflection profiles on image 269
       Modelled   210 /   214 reflection profiles on image 270
       Modelled   202 /   212 reflection profiles on image 271
       Modelled   194 /   205 reflection profiles on image 272
       Modelled   202 /   212 reflection profiles on image 273
       Modelled   226 /   233 reflection profiles on image 274
       Modelled   193 /   201 reflection profiles on image 275
       Modelled   221 /   225 reflection profiles on image 276
       Modelled   199 /   208 reflection profiles on image 277
       Modelled   220 /   236 reflection profiles on image 278
       Modelled   207 /   216 reflection profiles on image 279
       Modelled   186 /   202 reflection profiles on image 280
       Modelled   180 /   191 reflection profiles on image 281
       Modelled   225 /   231 reflection profiles on image 282
       Modelled   195 /   200 reflection profiles on image 283
       Modelled   204 /   212 reflection profiles on image 284
       Modelled    93 /    97 reflection profiles on image 285
       Modelled    36 /    37 reflection profiles on image 286
       Modelled    17 /    17 reflection profiles on image 287
       Modelled    14 /    15 reflection profiles on image 288
       Modelled     8 /     8 reflection profiles on image 289
       Modelled     4 /     4 reflection profiles on image 290
       Modelled     1 /     1 reflection profiles on image 292
       Modelled     3 /    14 reflection profiles on image 293
       Beginning modelling job 13

       Frames: 273 -> 315

       Number of reflections
        Partial:     37
        Full:        4512
        In ice ring: 0
        Total:       4549

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       273 [3   ]:
       274 [6   ]:
       275 [6   ]:
       276 [7   ]:
       277 [11  ]:
       278 [14  ]:
       279 [31  ]: *
       280 [56  ]: ***
       281 [91  ]: *****
       282 [292 ]: ****************
       283 [514 ]: *****************************
       284 [760 ]: *******************************************
       285 [981 ]: ********************************************************
       286 [1054]: ************************************************************
       287 [1087]: **************************************************************
       288 [1089]: **************************************************************
       289 [1079]: *************************************************************
       290 [1079]: *************************************************************
       291 [1081]: *************************************************************
       292 [1079]: *************************************************************
       293 [1100]: **************************************************************
       294 [1106]: ***************************************************************
       295 [1083]: **************************************************************
       296 [1076]: *************************************************************
       297 [1083]: **************************************************************
       298 [1107]: ***************************************************************
       299 [1146]: *****************************************************************
       300 [1170]: *******************************************************************
       301 [1120]: ****************************************************************
       302 [1044]: ***********************************************************
       303 [820 ]: **********************************************
       304 [597 ]: **********************************
       305 [381 ]: *********************
       306 [189 ]: **********
       307 [95  ]: *****
       308 [56  ]: ***
       309 [40  ]: **
       310 [31  ]: *
       311 [27  ]: *
       312 [20  ]: *
       313 [19  ]: *
       314 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0120292 GB

       Modelled   137 /   140 reflection profiles on image 285
       Modelled   179 /   185 reflection profiles on image 286
       Modelled   194 /   200 reflection profiles on image 287
       Modelled   223 /   229 reflection profiles on image 288
       Modelled   203 /   211 reflection profiles on image 289
       Modelled   182 /   188 reflection profiles on image 290
       Modelled   210 /   215 reflection profiles on image 291
       Modelled   210 /   218 reflection profiles on image 292
       Modelled   192 /   210 reflection profiles on image 293
       Modelled   194 /   206 reflection profiles on image 294
       Modelled   194 /   200 reflection profiles on image 295
       Modelled   213 /   222 reflection profiles on image 296
       Modelled   199 /   205 reflection profiles on image 297
       Modelled   187 /   194 reflection profiles on image 298
       Modelled   180 /   185 reflection profiles on image 299
       Modelled   231 /   241 reflection profiles on image 300
       Modelled   223 /   234 reflection profiles on image 301
       Modelled   240 /   246 reflection profiles on image 302
       Modelled   218 /   223 reflection profiles on image 303
       Modelled   205 /   216 reflection profiles on image 304
       Modelled   185 /   192 reflection profiles on image 305
       Modelled    92 /    94 reflection profiles on image 306
       Modelled    38 /    39 reflection profiles on image 307
       Modelled    15 /    16 reflection profiles on image 308
       Modelled     8 /     9 reflection profiles on image 309
       Modelled     4 /     4 reflection profiles on image 310
       Modelled     7 /     7 reflection profiles on image 311
       Modelled     1 /     1 reflection profiles on image 312
       Modelled     2 /     3 reflection profiles on image 313
       Modelled     0 /    16 reflection profiles on image 314
       Beginning modelling job 14

       Frames: 294 -> 336

       Number of reflections
        Partial:     37
        Full:        4462
        In ice ring: 0
        Total:       4499

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       295 [4   ]:
       296 [6   ]:
       297 [11  ]:
       298 [18  ]: *
       299 [28  ]: *
       300 [39  ]: **
       301 [55  ]: ***
       302 [108 ]: ******
       303 [291 ]: *****************
       304 [489 ]: ****************************
       305 [715 ]: ******************************************
       306 [923 ]: ******************************************************
       307 [1004]: ***********************************************************
       308 [1034]: *************************************************************
       309 [1030]: ************************************************************
       310 [1048]: **************************************************************
       311 [1080]: ***************************************************************
       312 [1121]: ******************************************************************
       313 [1118]: ******************************************************************
       314 [1104]: *****************************************************************
       315 [1108]: *****************************************************************
       316 [1127]: ******************************************************************
       317 [1132]: *******************************************************************
       318 [1096]: ****************************************************************
       319 [1118]: ******************************************************************
       320 [1102]: *****************************************************************
       321 [1091]: ****************************************************************
       322 [1058]: **************************************************************
       323 [1049]: **************************************************************
       324 [856 ]: **************************************************
       325 [633 ]: *************************************
       326 [417 ]: ************************
       327 [204 ]: ************
       328 [111 ]: ******
       329 [64  ]: ***
       330 [49  ]: **
       331 [40  ]: **
       332 [30  ]: *
       333 [23  ]: *
       334 [17  ]: *
       335 [14  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0123161 GB

       Modelled   120 /   123 reflection profiles on image 306
       Modelled   165 /   172 reflection profiles on image 307
       Modelled   201 /   210 reflection profiles on image 308
       Modelled   196 /   201 reflection profiles on image 309
       Modelled   194 /   197 reflection profiles on image 310
       Modelled   172 /   179 reflection profiles on image 311
       Modelled   198 /   206 reflection profiles on image 312
       Modelled   208 /   215 reflection profiles on image 313
       Modelled   206 /   221 reflection profiles on image 314
       Modelled   211 /   214 reflection profiles on image 315
       Modelled   194 /   202 reflection profiles on image 316
       Modelled   231 /   239 reflection profiles on image 317
       Modelled   183 /   190 reflection profiles on image 318
       Modelled   215 /   226 reflection profiles on image 319
       Modelled   201 /   207 reflection profiles on image 320
       Modelled   218 /   226 reflection profiles on image 321
       Modelled   192 /   197 reflection profiles on image 322
       Modelled   210 /   218 reflection profiles on image 323
       Modelled   215 /   223 reflection profiles on image 324
       Modelled   208 /   216 reflection profiles on image 325
       Modelled   203 /   213 reflection profiles on image 326
       Modelled    91 /    93 reflection profiles on image 327
       Modelled    44 /    47 reflection profiles on image 328
       Modelled    15 /    15 reflection profiles on image 329
       Modelled     9 /     9 reflection profiles on image 330
       Modelled     9 /    10 reflection profiles on image 331
       Modelled     5 /     7 reflection profiles on image 332
       Modelled     6 /     6 reflection profiles on image 333
       Modelled     3 /     3 reflection profiles on image 334
       Modelled     2 /    14 reflection profiles on image 335
       Beginning modelling job 15

       Frames: 315 -> 357

       Number of reflections
        Partial:     37
        Full:        4467
        In ice ring: 0
        Total:       4504

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       315 [2   ]:
       316 [2   ]:
       317 [2   ]:
       318 [7   ]:
       319 [9   ]:
       320 [14  ]:
       321 [29  ]: *
       322 [55  ]: ***
       323 [98  ]: *****
       324 [295 ]: *****************
       325 [510 ]: ******************************
       326 [731 ]: *******************************************
       327 [959 ]: *********************************************************
       328 [1016]: *************************************************************
       329 [1091]: *****************************************************************
       330 [1092]: *****************************************************************
       331 [1102]: ******************************************************************
       332 [1082]: *****************************************************************
       333 [1089]: *****************************************************************
       334 [1063]: ***************************************************************
       335 [1058]: ***************************************************************
       336 [1099]: ******************************************************************
       337 [1097]: *****************************************************************
       338 [1114]: *******************************************************************
       339 [1096]: *****************************************************************
       340 [1106]: ******************************************************************
       341 [1100]: ******************************************************************
       342 [1085]: *****************************************************************
       343 [1055]: ***************************************************************
       344 [1023]: *************************************************************
       345 [836 ]: **************************************************
       346 [609 ]: ************************************
       347 [381 ]: **********************
       348 [180 ]: **********
       349 [90  ]: *****
       350 [60  ]: ***
       351 [39  ]: **
       352 [34  ]: **
       353 [26  ]: *
       354 [23  ]: *
       355 [20  ]: *
       356 [15  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0117427 GB

       Modelled   134 /   138 reflection profiles on image 327
       Modelled   148 /   157 reflection profiles on image 328
       Modelled   198 /   206 reflection profiles on image 329
       Modelled   208 /   214 reflection profiles on image 330
       Modelled   206 /   212 reflection profiles on image 331
       Modelled   197 /   202 reflection profiles on image 332
       Modelled   218 /   224 reflection profiles on image 333
       Modelled   196 /   209 reflection profiles on image 334
       Modelled   195 /   214 reflection profiles on image 335
       Modelled   201 /   205 reflection profiles on image 336
       Modelled   190 /   205 reflection profiles on image 337
       Modelled   212 /   220 reflection profiles on image 338
       Modelled   197 /   205 reflection profiles on image 339
       Modelled   207 /   215 reflection profiles on image 340
       Modelled   202 /   207 reflection profiles on image 341
       Modelled   220 /   228 reflection profiles on image 342
       Modelled   190 /   197 reflection profiles on image 343
       Modelled   206 /   210 reflection profiles on image 344
       Modelled   221 /   227 reflection profiles on image 345
       Modelled   223 /   228 reflection profiles on image 346
       Modelled   197 /   201 reflection profiles on image 347
       Modelled    87 /    90 reflection profiles on image 348
       Modelled    28 /    30 reflection profiles on image 349
       Modelled    20 /    21 reflection profiles on image 350
       Modelled     4 /     5 reflection profiles on image 351
       Modelled     8 /     8 reflection profiles on image 352
       Modelled     3 /     3 reflection profiles on image 353
       Modelled     3 /     3 reflection profiles on image 354
       Modelled     4 /     5 reflection profiles on image 355
       Modelled     2 /    15 reflection profiles on image 356
       Beginning modelling job 16

       Frames: 336 -> 378

       Number of reflections
        Partial:     40
        Full:        4528
        In ice ring: 0
        Total:       4568

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       336 [3   ]:
       337 [7   ]:
       338 [9   ]:
       339 [16  ]:
       340 [19  ]: *
       341 [27  ]: *
       342 [44  ]: **
       343 [72  ]: ****
       344 [127 ]: *******
       345 [324 ]: ******************
       346 [548 ]: *******************************
       347 [774 ]: ********************************************
       348 [1008]: **********************************************************
       349 [1093]: ***************************************************************
       350 [1135]: *****************************************************************
       351 [1115]: ****************************************************************
       352 [1078]: **************************************************************
       353 [1082]: **************************************************************
       354 [1042]: ************************************************************
       355 [1069]: *************************************************************
       356 [1070]: *************************************************************
       357 [1102]: ***************************************************************
       358 [1101]: ***************************************************************
       359 [1089]: **************************************************************
       360 [1098]: ***************************************************************
       361 [1134]: *****************************************************************
       362 [1160]: *******************************************************************
       363 [1150]: ******************************************************************
       364 [1100]: ***************************************************************
       365 [1048]: ************************************************************
       366 [860 ]: *************************************************
       367 [644 ]: *************************************
       368 [419 ]: ************************
       369 [195 ]: ***********
       370 [96  ]: *****
       371 [64  ]: ***
       372 [47  ]: **
       373 [32  ]: *
       374 [25  ]: *
       375 [20  ]: *
       376 [17  ]:
       377 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0123561 GB

       Modelled   123 /   126 reflection profiles on image 348
       Modelled   164 /   169 reflection profiles on image 349
       Modelled   216 /   223 reflection profiles on image 350
       Modelled   228 /   234 reflection profiles on image 351
       Modelled   206 /   215 reflection profiles on image 352
       Modelled   211 /   218 reflection profiles on image 353
       Modelled   195 /   200 reflection profiles on image 354
       Modelled   195 /   204 reflection profiles on image 355
       Modelled   182 /   204 reflection profiles on image 356
       Modelled   194 /   198 reflection profiles on image 357
       Modelled   198 /   207 reflection profiles on image 358
       Modelled   215 /   219 reflection profiles on image 359
       Modelled   191 /   198 reflection profiles on image 360
       Modelled   188 /   190 reflection profiles on image 361
       Modelled   207 /   215 reflection profiles on image 362
       Modelled   248 /   263 reflection profiles on image 363
       Modelled   203 /   211 reflection profiles on image 364
       Modelled   207 /   214 reflection profiles on image 365
       Modelled   211 /   216 reflection profiles on image 366
       Modelled   217 /   225 reflection profiles on image 367
       Modelled   215 /   224 reflection profiles on image 368
       Modelled    94 /    99 reflection profiles on image 369
       Modelled    31 /    32 reflection profiles on image 370
       Modelled    15 /    17 reflection profiles on image 371
       Modelled    14 /    15 reflection profiles on image 372
       Modelled     7 /     7 reflection profiles on image 373
       Modelled     5 /     5 reflection profiles on image 374
       Modelled     3 /     3 reflection profiles on image 375
       Modelled     1 /     1 reflection profiles on image 376
       Modelled     2 /    16 reflection profiles on image 377
       Beginning modelling job 17

       Frames: 357 -> 399

       Number of reflections
        Partial:     38
        Full:        4519
        In ice ring: 0
        Total:       4557

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       357 [1   ]:
       358 [3   ]:
       359 [6   ]:
       360 [15  ]:
       361 [22  ]: *
       362 [35  ]: **
       363 [43  ]: **
       364 [66  ]: ***
       365 [96  ]: *****
       366 [319 ]: ******************
       367 [556 ]: ********************************
       368 [772 ]: *********************************************
       369 [946 ]: ********************************************************
       370 [1025]: ************************************************************
       371 [1021]: ************************************************************
       372 [1054]: **************************************************************
       373 [1055]: **************************************************************
       374 [1081]: ****************************************************************
       375 [1110]: *****************************************************************
       376 [1103]: *****************************************************************
       377 [1112]: *****************************************************************
       378 [1122]: ******************************************************************
       379 [1129]: *******************************************************************
       380 [1086]: ****************************************************************
       381 [1062]: ***************************************************************
       382 [1076]: ***************************************************************
       383 [1092]: ****************************************************************
       384 [1097]: *****************************************************************
       385 [1097]: *****************************************************************
       386 [1084]: ****************************************************************
       387 [881 ]: ****************************************************
       388 [663 ]: ***************************************
       389 [430 ]: *************************
       390 [194 ]: ***********
       391 [92  ]: *****
       392 [60  ]: ***
       393 [42  ]: **
       394 [28  ]: *
       395 [24  ]: *
       396 [20  ]: *
       397 [16  ]:
       398 [14  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0121727 GB

       Modelled   142 /   145 reflection profiles on image 369
       Modelled   186 /   191 reflection profiles on image 370
       Modelled   201 /   206 reflection profiles on image 371
       Modelled   190 /   197 reflection profiles on image 372
       Modelled   178 /   183 reflection profiles on image 373
       Modelled   189 /   194 reflection profiles on image 374
       Modelled   219 /   224 reflection profiles on image 375
       Modelled   201 /   205 reflection profiles on image 376
       Modelled   187 /   208 reflection profiles on image 377
       Modelled   219 /   224 reflection profiles on image 378
       Modelled   222 /   232 reflection profiles on image 379
       Modelled   204 /   217 reflection profiles on image 380
       Modelled   189 /   199 reflection profiles on image 381
       Modelled   190 /   198 reflection profiles on image 382
       Modelled   189 /   198 reflection profiles on image 383
       Modelled   205 /   211 reflection profiles on image 384
       Modelled   206 /   214 reflection profiles on image 385
       Modelled   224 /   230 reflection profiles on image 386
       Modelled   212 /   218 reflection profiles on image 387
       Modelled   225 /   233 reflection profiles on image 388
       Modelled   230 /   236 reflection profiles on image 389
       Modelled    97 /   102 reflection profiles on image 390
       Modelled    30 /    32 reflection profiles on image 391
       Modelled    18 /    18 reflection profiles on image 392
       Modelled    13 /    14 reflection profiles on image 393
       Modelled     4 /     4 reflection profiles on image 394
       Modelled     4 /     4 reflection profiles on image 395
       Modelled     3 /     4 reflection profiles on image 396
       Modelled     2 /     2 reflection profiles on image 397
       Modelled     0 /    14 reflection profiles on image 398
       Beginning modelling job 18

       Frames: 378 -> 420

       Number of reflections
        Partial:     40
        Full:        4578
        In ice ring: 0
        Total:       4618

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       378 [2   ]:
       379 [5   ]:
       380 [6   ]:
       381 [9   ]:
       382 [16  ]:
       383 [24  ]: *
       384 [45  ]: **
       385 [66  ]: ***
       386 [122 ]: *******
       387 [341 ]: *******************
       388 [599 ]: **********************************
       389 [827 ]: ***********************************************
       390 [1050]: ************************************************************
       391 [1155]: *******************************************************************
       392 [1108]: ****************************************************************
       393 [1112]: ****************************************************************
       394 [1111]: ****************************************************************
       395 [1081]: **************************************************************
       396 [1124]: *****************************************************************
       397 [1126]: *****************************************************************
       398 [1115]: ****************************************************************
       399 [1124]: *****************************************************************
       400 [1106]: ****************************************************************
       401 [1102]: ***************************************************************
       402 [1097]: ***************************************************************
       403 [1105]: ****************************************************************
       404 [1094]: ***************************************************************
       405 [1090]: ***************************************************************
       406 [1091]: ***************************************************************
       407 [1061]: *************************************************************
       408 [854 ]: *************************************************
       409 [626 ]: ************************************
       410 [430 ]: ************************
       411 [198 ]: ***********
       412 [105 ]: ******
       413 [61  ]: ***
       414 [46  ]: **
       415 [35  ]: **
       416 [32  ]: *
       417 [26  ]: *
       418 [20  ]: *
       419 [19  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0123498 GB

       Modelled   127 /   130 reflection profiles on image 390
       Modelled   221 /   229 reflection profiles on image 391
       Modelled   198 /   202 reflection profiles on image 392
       Modelled   209 /   217 reflection profiles on image 393
       Modelled   225 /   233 reflection profiles on image 394
       Modelled   173 /   181 reflection profiles on image 395
       Modelled   217 /   222 reflection profiles on image 396
       Modelled   212 /   219 reflection profiles on image 397
       Modelled   206 /   225 reflection profiles on image 398
       Modelled   216 /   222 reflection profiles on image 399
       Modelled   196 /   207 reflection profiles on image 400
       Modelled   210 /   217 reflection profiles on image 401
       Modelled   194 /   201 reflection profiles on image 402
       Modelled   204 /   214 reflection profiles on image 403
       Modelled   198 /   202 reflection profiles on image 404
       Modelled   197 /   205 reflection profiles on image 405
       Modelled   200 /   209 reflection profiles on image 406
       Modelled   221 /   229 reflection profiles on image 407
       Modelled   216 /   228 reflection profiles on image 408
       Modelled   188 /   196 reflection profiles on image 409
       Modelled   221 /   232 reflection profiles on image 410
       Modelled    93 /    93 reflection profiles on image 411
       Modelled    42 /    44 reflection profiles on image 412
       Modelled    15 /    15 reflection profiles on image 413
       Modelled     9 /    11 reflection profiles on image 414
       Modelled     3 /     3 reflection profiles on image 415
       Modelled     3 /     6 reflection profiles on image 416
       Modelled     6 /     6 reflection profiles on image 417
       Modelled     1 /     1 reflection profiles on image 418
       Modelled     5 /    19 reflection profiles on image 419
       Beginning modelling job 19

       Frames: 399 -> 441

       Number of reflections
        Partial:     41
        Full:        4488
        In ice ring: 0
        Total:       4529

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       399 [2   ]:
       400 [3   ]:
       401 [4   ]:
       402 [8   ]:
       403 [15  ]:
       404 [26  ]: *
       405 [33  ]: *
       406 [58  ]: ***
       407 [106 ]: *****
       408 [299 ]: ****************
       409 [520 ]: *****************************
       410 [769 ]: *******************************************
       411 [961 ]: *****************************************************
       412 [1040]: **********************************************************
       413 [1083]: ************************************************************
       414 [1103]: *************************************************************
       415 [1146]: ****************************************************************
       416 [1145]: ****************************************************************
       417 [1154]: ****************************************************************
       418 [1181]: ******************************************************************
       419 [1174]: *****************************************************************
       420 [1195]: *******************************************************************
       421 [1159]: ****************************************************************
       422 [1133]: ***************************************************************
       423 [1097]: *************************************************************
       424 [1087]: ************************************************************
       425 [1067]: ***********************************************************
       426 [1031]: *********************************************************
       427 [1015]: ********************************************************
       428 [959 ]: *****************************************************
       429 [789 ]: ********************************************
       430 [578 ]: ********************************
       431 [387 ]: *********************
       432 [154 ]: ********
       433 [89  ]: ****
       434 [60  ]: ***
       435 [39  ]: **
       436 [33  ]: *
       437 [27  ]: *
       438 [21  ]: *
       439 [20  ]: *
       440 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0129768 GB

       Modelled   121 /   125 reflection profiles on image 411
       Modelled   169 /   175 reflection profiles on image 412
       Modelled   214 /   222 reflection profiles on image 413
       Modelled   177 /   184 reflection profiles on image 414
       Modelled   193 /   199 reflection profiles on image 415
       Modelled   219 /   229 reflection profiles on image 416
       Modelled   190 /   198 reflection profiles on image 417
       Modelled   218 /   221 reflection profiles on image 418
       Modelled   219 /   236 reflection profiles on image 419
       Modelled   214 /   222 reflection profiles on image 420
       Modelled   232 /   238 reflection profiles on image 421
       Modelled   218 /   227 reflection profiles on image 422
       Modelled   219 /   231 reflection profiles on image 423
       Modelled   210 /   218 reflection profiles on image 424
       Modelled   212 /   223 reflection profiles on image 425
       Modelled   176 /   185 reflection profiles on image 426
       Modelled   196 /   207 reflection profiles on image 427
       Modelled   192 /   200 reflection profiles on image 428
       Modelled   209 /   211 reflection profiles on image 429
       Modelled   182 /   191 reflection profiles on image 430
       Modelled   228 /   233 reflection profiles on image 431
       Modelled    64 /    65 reflection profiles on image 432
       Modelled    25 /    29 reflection profiles on image 433
       Modelled    20 /    21 reflection profiles on image 434
       Modelled     6 /     6 reflection profiles on image 435
       Modelled     6 /     6 reflection profiles on image 436
       Modelled     5 /     6 reflection profiles on image 437
       Modelled     1 /     1 reflection profiles on image 438
       Modelled     3 /     4 reflection profiles on image 439
       Modelled     3 /    16 reflection profiles on image 440
       Beginning modelling job 20

       Frames: 420 -> 462

       Number of reflections
        Partial:     40
        Full:        4431
        In ice ring: 0
        Total:       4471

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       420 [1   ]:
       421 [6   ]:
       422 [9   ]:
       423 [13  ]:
       424 [15  ]:
       425 [24  ]: *
       426 [36  ]: **
       427 [57  ]: ***
       428 [107 ]: ******
       429 [306 ]: ******************
       430 [510 ]: ******************************
       431 [718 ]: ******************************************
       432 [948 ]: ********************************************************
       433 [1022]: *************************************************************
       434 [1055]: **************************************************************
       435 [1072]: ****************************************************************
       436 [1073]: ****************************************************************
       437 [1049]: **************************************************************
       438 [1052]: **************************************************************
       439 [1080]: ****************************************************************
       440 [1085]: ****************************************************************
       441 [1120]: ******************************************************************
       442 [1122]: *******************************************************************
       443 [1110]: ******************************************************************
       444 [1090]: *****************************************************************
       445 [1091]: *****************************************************************
       446 [1083]: ****************************************************************
       447 [1067]: ***************************************************************
       448 [1028]: *************************************************************
       449 [1019]: ************************************************************
       450 [838 ]: **************************************************
       451 [624 ]: *************************************
       452 [429 ]: *************************
       453 [197 ]: ***********
       454 [95  ]: *****
       455 [61  ]: ***
       456 [49  ]: **
       457 [32  ]: *
       458 [24  ]: *
       459 [19  ]: *
       460 [19  ]: *
       461 [16  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.012184 GB

       Modelled   136 /   138 reflection profiles on image 432
       Modelled   165 /   169 reflection profiles on image 433
       Modelled   198 /   201 reflection profiles on image 434
       Modelled   175 /   185 reflection profiles on image 435
       Modelled   222 /   230 reflection profiles on image 436
       Modelled   199 /   203 reflection profiles on image 437
       Modelled   192 /   200 reflection profiles on image 438
       Modelled   193 /   196 reflection profiles on image 439
       Modelled   204 /   220 reflection profiles on image 440
       Modelled   196 /   202 reflection profiles on image 441
       Modelled   211 /   219 reflection profiles on image 442
       Modelled   221 /   223 reflection profiles on image 443
       Modelled   201 /   209 reflection profiles on image 444
       Modelled   201 /   208 reflection profiles on image 445
       Modelled   198 /   204 reflection profiles on image 446
       Modelled   213 /   221 reflection profiles on image 447
       Modelled   185 /   195 reflection profiles on image 448
       Modelled   202 /   210 reflection profiles on image 449
       Modelled   213 /   214 reflection profiles on image 450
       Modelled   190 /   195 reflection profiles on image 451
       Modelled   226 /   232 reflection profiles on image 452
       Modelled    95 /   102 reflection profiles on image 453
       Modelled    33 /    34 reflection profiles on image 454
       Modelled    12 /    12 reflection profiles on image 455
       Modelled    17 /    17 reflection profiles on image 456
       Modelled     6 /     8 reflection profiles on image 457
       Modelled     5 /     5 reflection profiles on image 458
       Modelled     3 /     3 reflection profiles on image 460
       Modelled     3 /    16 reflection profiles on image 461
       Beginning modelling job 21

       Frames: 441 -> 483

       Number of reflections
        Partial:     32
        Full:        4479
        In ice ring: 0
        Total:       4511

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       442 [6   ]:
       443 [11  ]:
       444 [15  ]:
       445 [18  ]: *
       446 [30  ]: *
       447 [39  ]: **
       448 [65  ]: ***
       449 [116 ]: ******
       450 [304 ]: *****************
       451 [517 ]: *****************************
       452 [707 ]: ****************************************
       453 [947 ]: *****************************************************
       454 [1008]: *********************************************************
       455 [1041]: ***********************************************************
       456 [1057]: ************************************************************
       457 [1076]: *************************************************************
       458 [1098]: **************************************************************
       459 [1129]: ****************************************************************
       460 [1147]: *****************************************************************
       461 [1170]: ******************************************************************
       462 [1178]: *******************************************************************
       463 [1171]: ******************************************************************
       464 [1168]: ******************************************************************
       465 [1148]: *****************************************************************
       466 [1100]: **************************************************************
       467 [1072]: ************************************************************
       468 [1056]: ************************************************************
       469 [1060]: ************************************************************
       470 [1016]: *********************************************************
       471 [803 ]: *********************************************
       472 [585 ]: *********************************
       473 [392 ]: **********************
       474 [178 ]: **********
       475 [94  ]: *****
       476 [62  ]: ***
       477 [46  ]: **
       478 [34  ]: *
       479 [24  ]: *
       480 [17  ]:
       481 [15  ]:
       482 [13  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0126419 GB

       Modelled   123 /   127 reflection profiles on image 453
       Modelled   172 /   182 reflection profiles on image 454
       Modelled   189 /   199 reflection profiles on image 455
       Modelled   181 /   184 reflection profiles on image 456
       Modelled   182 /   187 reflection profiles on image 457
       Modelled   195 /   200 reflection profiles on image 458
       Modelled   205 /   210 reflection profiles on image 459
       Modelled   201 /   207 reflection profiles on image 460
       Modelled   214 /   233 reflection profiles on image 461
       Modelled   202 /   217 reflection profiles on image 462
       Modelled   215 /   221 reflection profiles on image 463
       Modelled   219 /   227 reflection profiles on image 464
       Modelled   238 /   243 reflection profiles on image 465
       Modelled   232 /   240 reflection profiles on image 466
       Modelled   192 /   202 reflection profiles on image 467
       Modelled   185 /   192 reflection profiles on image 468
       Modelled   202 /   206 reflection profiles on image 469
       Modelled   227 /   231 reflection profiles on image 470
       Modelled   210 /   218 reflection profiles on image 471
       Modelled   185 /   193 reflection profiles on image 472
       Modelled   208 /   214 reflection profiles on image 473
       Modelled    79 /    84 reflection profiles on image 474
       Modelled    32 /    32 reflection profiles on image 475
       Modelled    14 /    16 reflection profiles on image 476
       Modelled    10 /    12 reflection profiles on image 477
       Modelled     9 /    10 reflection profiles on image 478
       Modelled     7 /     7 reflection profiles on image 479
       Modelled     2 /     2 reflection profiles on image 480
       Modelled     2 /     2 reflection profiles on image 481
       Modelled     3 /    13 reflection profiles on image 482
       Beginning modelling job 22

       Frames: 462 -> 504

       Number of reflections
        Partial:     76
        Full:        4367
        In ice ring: 0
        Total:       4443

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       462 [1   ]:
       463 [7   ]:
       464 [11  ]:
       465 [14  ]:
       466 [19  ]: *
       467 [23  ]: *
       468 [44  ]: **
       469 [71  ]: ****
       470 [114 ]: ******
       471 [275 ]: ****************
       472 [505 ]: *****************************
       473 [719 ]: *****************************************
       474 [932 ]: ******************************************************
       475 [1046]: ************************************************************
       476 [1082]: ***************************************************************
       477 [1094]: ***************************************************************
       478 [1086]: ***************************************************************
       479 [1094]: ***************************************************************
       480 [1059]: *************************************************************
       481 [1106]: ****************************************************************
       482 [1100]: ****************************************************************
       483 [1119]: *****************************************************************
       484 [1143]: ******************************************************************
       485 [1150]: *******************************************************************
       486 [1142]: ******************************************************************
       487 [1127]: *****************************************************************
       488 [1124]: *****************************************************************
       489 [1077]: **************************************************************
       490 [1052]: *************************************************************
       491 [994 ]: *********************************************************
       492 [800 ]: **********************************************
       493 [602 ]: ***********************************
       494 [390 ]: **********************
       495 [188 ]: **********
       496 [111 ]: ******
       497 [88  ]: *****
       498 [61  ]: ***
       499 [48  ]: **
       500 [45  ]: **
       501 [42  ]: **
       502 [37  ]: **
       503 [35  ]: **

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0125611 GB

       Modelled   101 /   103 reflection profiles on image 474
       Modelled   165 /   170 reflection profiles on image 475
       Modelled   180 /   185 reflection profiles on image 476
       Modelled   197 /   203 reflection profiles on image 477
       Modelled   202 /   212 reflection profiles on image 478
       Modelled   219 /   224 reflection profiles on image 479
       Modelled   164 /   171 reflection profiles on image 480
       Modelled   203 /   211 reflection profiles on image 481
       Modelled   204 /   241 reflection profiles on image 482
       Modelled   184 /   191 reflection profiles on image 483
       Modelled   203 /   208 reflection profiles on image 484
       Modelled   227 /   232 reflection profiles on image 485
       Modelled   201 /   208 reflection profiles on image 486
       Modelled   193 /   201 reflection profiles on image 487
       Modelled   226 /   235 reflection profiles on image 488
       Modelled   212 /   217 reflection profiles on image 489
       Modelled   209 /   215 reflection profiles on image 490
       Modelled   208 /   216 reflection profiles on image 491
       Modelled   192 /   198 reflection profiles on image 492
       Modelled   202 /   212 reflection profiles on image 493
       Modelled   195 /   202 reflection profiles on image 494
       Modelled    75 /    77 reflection profiles on image 495
       Modelled    23 /    23 reflection profiles on image 496
       Modelled    26 /    27 reflection profiles on image 497
       Modelled    11 /    13 reflection profiles on image 498
       Modelled     3 /     3 reflection profiles on image 499
       Modelled     2 /     3 reflection profiles on image 500
       Modelled     3 /     5 reflection profiles on image 501
       Modelled     1 /     2 reflection profiles on image 502
       Modelled     5 /    35 reflection profiles on image 503
       Beginning modelling job 23

       Frames: 483 -> 525

       Number of reflections
        Partial:     46
        Full:        4150
        In ice ring: 0
        Total:       4196

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       483 [1   ]:
       484 [5   ]:
       485 [7   ]:
       486 [13  ]:
       487 [18  ]: *
       488 [32  ]: *
       489 [45  ]: **
       490 [63  ]: ***
       491 [122 ]: *******
       492 [311 ]: *****************
       493 [544 ]: *******************************
       494 [754 ]: *******************************************
       495 [958 ]: *******************************************************
       496 [1031]: ***********************************************************
       497 [1065]: *************************************************************
       498 [1087]: **************************************************************
       499 [1099]: ***************************************************************
       500 [1102]: ***************************************************************
       501 [1100]: ***************************************************************
       502 [1115]: ****************************************************************
       503 [1080]: **************************************************************
       504 [1121]: ****************************************************************
       505 [1141]: *****************************************************************
       506 [1165]: *******************************************************************
       507 [1153]: ******************************************************************
       508 [1106]: ***************************************************************
       509 [1107]: ***************************************************************
       510 [1055]: ************************************************************
       511 [944 ]: ******************************************************
       512 [752 ]: *******************************************
       513 [512 ]: *****************************
       514 [295 ]: ****************
       515 [133 ]: *******
       516 [71  ]: ****
       517 [46  ]: **
       518 [33  ]: *
       519 [19  ]: *
       520 [13  ]:
       521 [7   ]:
       522 [3   ]:
       523 [2   ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.012759 GB

       Modelled   122 /   128 reflection profiles on image 495
       Modelled   183 /   187 reflection profiles on image 496
       Modelled   196 /   204 reflection profiles on image 497
       Modelled   210 /   217 reflection profiles on image 498
       Modelled   190 /   191 reflection profiles on image 499
       Modelled   183 /   192 reflection profiles on image 500
       Modelled   194 /   202 reflection profiles on image 501
       Modelled   229 /   237 reflection profiles on image 502
       Modelled   184 /   207 reflection profiles on image 503
       Modelled   188 /   201 reflection profiles on image 504
       Modelled   196 /   203 reflection profiles on image 505
       Modelled   195 /   202 reflection profiles on image 506
       Modelled   210 /   220 reflection profiles on image 507
       Modelled   211 /   217 reflection profiles on image 508
       Modelled   227 /   232 reflection profiles on image 509
       Modelled   202 /   212 reflection profiles on image 510
       Modelled   187 /   192 reflection profiles on image 511
       Modelled   229 /   240 reflection profiles on image 512
       Modelled   203 /   217 reflection profiles on image 513
       Modelled   157 /   162 reflection profiles on image 514
       Modelled    56 /    62 reflection profiles on image 515
       Modelled    21 /    25 reflection profiles on image 516
       Modelled    11 /    13 reflection profiles on image 517
       Modelled    11 /    14 reflection profiles on image 518
       Modelled     5 /     6 reflection profiles on image 519
       Modelled     6 /     6 reflection profiles on image 520
       Modelled     3 /     4 reflection profiles on image 521
       Modelled     1 /     1 reflection profiles on image 522
       Modelled     2 /     2 reflection profiles on image 523
       Beginning modelling job 24

       Frames: 504 -> 540

       Number of reflections
        Partial:     434
        Full:        5246
        In ice ring: 0
        Total:       5680

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       504 [18  ]: *
       505 [20  ]: *
       506 [23  ]: *
       507 [27  ]: *
       508 [38  ]: **
       509 [56  ]: ***
       510 [84  ]: ****
       511 [182 ]: **********
       512 [386 ]: **********************
       513 [608 ]: **********************************
       514 [826 ]: ***********************************************
       515 [1029]: **********************************************************
       516 [1069]: *************************************************************
       517 [1072]: *************************************************************
       518 [1085]: **************************************************************
       519 [1101]: **************************************************************
       520 [1124]: ****************************************************************
       521 [1125]: ****************************************************************
       522 [1126]: ****************************************************************
       523 [1142]: *****************************************************************
       524 [1162]: ******************************************************************
       525 [1172]: *******************************************************************
       526 [1154]: *****************************************************************
       527 [1133]: ****************************************************************
       528 [1132]: ****************************************************************
       529 [1116]: ***************************************************************
       530 [1138]: *****************************************************************
       531 [1143]: *****************************************************************
       532 [1123]: ****************************************************************
       533 [1075]: *************************************************************
       534 [1071]: *************************************************************
       535 [1048]: ***********************************************************
       536 [1034]: ***********************************************************
       537 [972 ]: *******************************************************
       538 [835 ]: ***********************************************
       539 [616 ]: ***********************************

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0127654 GB

       Modelled    22 /    22 reflection profiles on image 514
       Modelled   156 /   161 reflection profiles on image 515
       Modelled   188 /   198 reflection profiles on image 516
       Modelled   198 /   206 reflection profiles on image 517
       Modelled   194 /   203 reflection profiles on image 518
       Modelled   184 /   187 reflection profiles on image 519
       Modelled   212 /   220 reflection profiles on image 520
       Modelled   213 /   221 reflection profiles on image 521
       Modelled   187 /   191 reflection profiles on image 522
       Modelled   200 /   214 reflection profiles on image 523
       Modelled   194 /   216 reflection profiles on image 524
       Modelled   213 /   225 reflection profiles on image 525
       Modelled   219 /   228 reflection profiles on image 526
       Modelled   204 /   215 reflection profiles on image 527
       Modelled   223 /   228 reflection profiles on image 528
       Modelled   208 /   212 reflection profiles on image 529
       Modelled   198 /   203 reflection profiles on image 530
       Modelled   205 /   211 reflection profiles on image 531
       Modelled   236 /   248 reflection profiles on image 532
       Modelled   194 /   197 reflection profiles on image 533
       Modelled   215 /   225 reflection profiles on image 534
       Modelled   196 /   202 reflection profiles on image 535
       Modelled   203 /   209 reflection profiles on image 536
       Modelled   196 /   203 reflection profiles on image 537
       Modelled   216 /   219 reflection profiles on image 538
       Modelled   272 /   616 reflection profiles on image 539

       Summary of profile model
       -------------------------------------------------------------------
       ID | Profile | Created | X (px)  | Y (px)  | Z (im) | # reflections
       -------------------------------------------------------------------
       0  | 0       | True    | 410.50  | 421.17  | 15.88  | 7309
       0  | 1       | True    | 1231.50 | 421.17  | 15.88  | 9195
       0  | 2       | True    | 2052.50 | 421.17  | 15.88  | 7356
       0  | 3       | True    | 410.50  | 1263.50 | 15.88  | 10215
       0  | 4       | True    | 1231.50 | 1263.50 | 15.88  | 12524
       0  | 5       | True    | 2052.50 | 1263.50 | 15.88  | 10002
       0  | 6       | True    | 410.50  | 2105.83 | 15.88  | 6553
       0  | 7       | True    | 1231.50 | 2105.83 | 15.88  | 7981
       0  | 8       | True    | 2052.50 | 2105.83 | 15.88  | 6206
       0  | 9       | True    | 410.50  | 421.17  | 47.65  | 11065
       0  | 10      | True    | 1231.50 | 421.17  | 47.65  | 13951
       0  | 11      | True    | 2052.50 | 421.17  | 47.65  | 11094
       0  | 12      | True    | 410.50  | 1263.50 | 47.65  | 15431
       0  | 13      | True    | 1231.50 | 1263.50 | 47.65  | 18917
       0  | 14      | True    | 2052.50 | 1263.50 | 47.65  | 14981
       0  | 15      | True    | 410.50  | 2105.83 | 47.65  | 9928
       0  | 16      | True    | 1231.50 | 2105.83 | 47.65  | 12059
       0  | 17      | True    | 2052.50 | 2105.83 | 47.65  | 9282
       0  | 18      | True    | 410.50  | 421.17  | 79.41  | 11221
       0  | 19      | True    | 1231.50 | 421.17  | 79.41  | 14277
       0  | 20      | True    | 2052.50 | 421.17  | 79.41  | 11429
       0  | 21      | True    | 410.50  | 1263.50 | 79.41  | 15496
       0  | 22      | True    | 1231.50 | 1263.50 | 79.41  | 19141
       0  | 23      | True    | 2052.50 | 1263.50 | 79.41  | 15258
       0  | 24      | True    | 410.50  | 2105.83 | 79.41  | 9925
       0  | 25      | True    | 1231.50 | 2105.83 | 79.41  | 12131
       0  | 26      | True    | 2052.50 | 2105.83 | 79.41  | 9362
       0  | 27      | True    | 410.50  | 421.17  | 111.18 | 11324
       0  | 28      | True    | 1231.50 | 421.17  | 111.18 | 14384
       0  | 29      | True    | 2052.50 | 421.17  | 111.18 | 11493
       0  | 30      | True    | 410.50  | 1263.50 | 111.18 | 15596
       0  | 31      | True    | 1231.50 | 1263.50 | 111.18 | 19239
       0  | 32      | True    | 2052.50 | 1263.50 | 111.18 | 15276
       0  | 33      | True    | 410.50  | 2105.83 | 111.18 | 9996
       0  | 34      | True    | 1231.50 | 2105.83 | 111.18 | 12192
       0  | 35      | True    | 2052.50 | 2105.83 | 111.18 | 9359
       0  | 36      | True    | 410.50  | 421.17  | 142.94 | 11303
       0  | 37      | True    | 1231.50 | 421.17  | 142.94 | 14403
       0  | 38      | True    | 2052.50 | 421.17  | 142.94 | 11562
       0  | 39      | True    | 410.50  | 1263.50 | 142.94 | 15584
       0  | 40      | True    | 1231.50 | 1263.50 | 142.94 | 19309
       0  | 41      | True    | 2052.50 | 1263.50 | 142.94 | 15410
       0  | 42      | True    | 410.50  | 2105.83 | 142.94 | 9920
       0  | 43      | True    | 1231.50 | 2105.83 | 142.94 | 12180
       0  | 44      | True    | 2052.50 | 2105.83 | 142.94 | 9371
       0  | 45      | True    | 410.50  | 421.17  | 174.71 | 11299
       0  | 46      | True    | 1231.50 | 421.17  | 174.71 | 14355
       0  | 47      | True    | 2052.50 | 421.17  | 174.71 | 11518
       0  | 48      | True    | 410.50  | 1263.50 | 174.71 | 15669
       0  | 49      | True    | 1231.50 | 1263.50 | 174.71 | 19360
       0  | 50      | True    | 2052.50 | 1263.50 | 174.71 | 15423
       0  | 51      | True    | 410.50  | 2105.83 | 174.71 | 10021
       0  | 52      | True    | 1231.50 | 2105.83 | 174.71 | 12277
       0  | 53      | True    | 2052.50 | 2105.83 | 174.71 | 9428
       0  | 54      | True    | 410.50  | 421.17  | 206.47 | 11205
       0  | 55      | True    | 1231.50 | 421.17  | 206.47 | 14264
       0  | 56      | True    | 2052.50 | 421.17  | 206.47 | 11464
       0  | 57      | True    | 410.50  | 1263.50 | 206.47 | 15619
       0  | 58      | True    | 1231.50 | 1263.50 | 206.47 | 19345
       0  | 59      | True    | 2052.50 | 1263.50 | 206.47 | 15465
       0  | 60      | True    | 410.50  | 2105.83 | 206.47 | 10017
       0  | 61      | True    | 1231.50 | 2105.83 | 206.47 | 12302
       0  | 62      | True    | 2052.50 | 2105.83 | 206.47 | 9477
       0  | 63      | True    | 410.50  | 421.17  | 238.24 | 11176
       0  | 64      | True    | 1231.50 | 421.17  | 238.24 | 14233
       0  | 65      | True    | 2052.50 | 421.17  | 238.24 | 11476
       0  | 66      | True    | 410.50  | 1263.50 | 238.24 | 15590
       0  | 67      | True    | 1231.50 | 1263.50 | 238.24 | 19332
       0  | 68      | True    | 2052.50 | 1263.50 | 238.24 | 15550
       0  | 69      | True    | 410.50  | 2105.83 | 238.24 | 10011
       0  | 70      | True    | 1231.50 | 2105.83 | 238.24 | 12305
       0  | 71      | True    | 2052.50 | 2105.83 | 238.24 | 9564
       0  | 72      | True    | 410.50  | 421.17  | 270.00 | 11087
       0  | 73      | True    | 1231.50 | 421.17  | 270.00 | 14190
       0  | 74      | True    | 2052.50 | 421.17  | 270.00 | 11421
       0  | 75      | True    | 410.50  | 1263.50 | 270.00 | 15567
       0  | 76      | True    | 1231.50 | 1263.50 | 270.00 | 19403
       0  | 77      | True    | 2052.50 | 1263.50 | 270.00 | 15585
       0  | 78      | True    | 410.50  | 2105.83 | 270.00 | 10055
       0  | 79      | True    | 1231.50 | 2105.83 | 270.00 | 12411
       0  | 80      | True    | 2052.50 | 2105.83 | 270.00 | 9638
       0  | 81      | True    | 410.50  | 421.17  | 301.76 | 11084
       0  | 82      | True    | 1231.50 | 421.17  | 301.76 | 14311
       0  | 83      | True    | 2052.50 | 421.17  | 301.76 | 11532
       0  | 84      | True    | 410.50  | 1263.50 | 301.76 | 15593
       0  | 85      | True    | 1231.50 | 1263.50 | 301.76 | 19578
       0  | 86      | True    | 2052.50 | 1263.50 | 301.76 | 15762
       0  | 87      | True    | 410.50  | 2105.83 | 301.76 | 10041
       0  | 88      | True    | 1231.50 | 2105.83 | 301.76 | 12478
       0  | 89      | True    | 2052.50 | 2105.83 | 301.76 | 9727
       0  | 90      | True    | 410.50  | 421.17  | 333.53 | 11185
       0  | 91      | True    | 1231.50 | 421.17  | 333.53 | 14452
       0  | 92      | True    | 2052.50 | 421.17  | 333.53 | 11625
       0  | 93      | True    | 410.50  | 1263.50 | 333.53 | 15680
       0  | 94      | True    | 1231.50 | 1263.50 | 333.53 | 19708
       0  | 95      | True    | 2052.50 | 1263.50 | 333.53 | 15805
       0  | 96      | True    | 410.50  | 2105.83 | 333.53 | 10083
       0  | 97      | True    | 1231.50 | 2105.83 | 333.53 | 12541
       0  | 98      | True    | 2052.50 | 2105.83 | 333.53 | 9729
       0  | 99      | True    | 410.50  | 421.17  | 365.29 | 11260
       0  | 100     | True    | 1231.50 | 421.17  | 365.29 | 14635
       0  | 101     | True    | 2052.50 | 421.17  | 365.29 | 11794
       0  | 102     | True    | 410.50  | 1263.50 | 365.29 | 15730
       0  | 103     | True    | 1231.50 | 1263.50 | 365.29 | 19898
       0  | 104     | True    | 2052.50 | 1263.50 | 365.29 | 15976
       0  | 105     | True    | 410.50  | 2105.83 | 365.29 | 9988
       0  | 106     | True    | 1231.50 | 2105.83 | 365.29 | 12492
       0  | 107     | True    | 2052.50 | 2105.83 | 365.29 | 9700
       0  | 108     | True    | 410.50  | 421.17  | 397.06 | 11284
       0  | 109     | True    | 1231.50 | 421.17  | 397.06 | 14634
       0  | 110     | True    | 2052.50 | 421.17  | 397.06 | 11718
       0  | 111     | True    | 410.50  | 1263.50 | 397.06 | 15662
       0  | 112     | True    | 1231.50 | 1263.50 | 397.06 | 19827
       0  | 113     | True    | 2052.50 | 1263.50 | 397.06 | 15810
       0  | 114     | True    | 410.50  | 2105.83 | 397.06 | 9874
       0  | 115     | True    | 1231.50 | 2105.83 | 397.06 | 12389
       0  | 116     | True    | 2052.50 | 2105.83 | 397.06 | 9570
       0  | 117     | True    | 410.50  | 421.17  | 428.82 | 11216
       0  | 118     | True    | 1231.50 | 421.17  | 428.82 | 14592
       0  | 119     | True    | 2052.50 | 421.17  | 428.82 | 11664
       0  | 120     | True    | 410.50  | 1263.50 | 428.82 | 15515
       0  | 121     | True    | 1231.50 | 1263.50 | 428.82 | 19753
       0  | 122     | True    | 2052.50 | 1263.50 | 428.82 | 15773
       0  | 123     | True    | 410.50  | 2105.83 | 428.82 | 9727
       0  | 124     | True    | 1231.50 | 2105.83 | 428.82 | 12294
       0  | 125     | True    | 2052.50 | 2105.83 | 428.82 | 9534
       0  | 126     | True    | 410.50  | 421.17  | 460.59 | 11193
       0  | 127     | True    | 1231.50 | 421.17  | 460.59 | 14441
       0  | 128     | True    | 2052.50 | 421.17  | 460.59 | 11495
       0  | 129     | True    | 410.50  | 1263.50 | 460.59 | 15427
       0  | 130     | True    | 1231.50 | 1263.50 | 460.59 | 19540
       0  | 131     | True    | 2052.50 | 1263.50 | 460.59 | 15558
       0  | 132     | True    | 410.50  | 2105.83 | 460.59 | 9642
       0  | 133     | True    | 1231.50 | 2105.83 | 460.59 | 12181
       0  | 134     | True    | 2052.50 | 2105.83 | 460.59 | 9419
       0  | 135     | True    | 410.50  | 421.17  | 492.35 | 10942
       0  | 136     | True    | 1231.50 | 421.17  | 492.35 | 14129
       0  | 137     | True    | 2052.50 | 421.17  | 492.35 | 11299
       0  | 138     | True    | 410.50  | 1263.50 | 492.35 | 15116
       0  | 139     | True    | 1231.50 | 1263.50 | 492.35 | 19163
       0  | 140     | True    | 2052.50 | 1263.50 | 492.35 | 15342
       0  | 141     | True    | 410.50  | 2105.83 | 492.35 | 9451
       0  | 142     | True    | 1231.50 | 2105.83 | 492.35 | 11941
       0  | 143     | True    | 2052.50 | 2105.83 | 492.35 | 9287
       0  | 144     | True    | 410.50  | 421.17  | 524.12 | 7203
       0  | 145     | True    | 1231.50 | 421.17  | 524.12 | 9286
       0  | 146     | True    | 2052.50 | 421.17  | 524.12 | 7441
       0  | 147     | True    | 410.50  | 1263.50 | 524.12 | 10001
       0  | 148     | True    | 1231.50 | 1263.50 | 524.12 | 12657
       0  | 149     | True    | 2052.50 | 1263.50 | 524.12 | 10138
       0  | 150     | True    | 410.50  | 2105.83 | 524.12 | 6263
       0  | 151     | True    | 1231.50 | 2105.83 | 524.12 | 7904
       0  | 152     | True    | 2052.50 | 2105.83 | 524.12 | 6146
       -------------------------------------------------------------------


       ----------------------------------
               Read time | 124.13 seconds
            Extract time |   1.07 seconds
        Pre-process time |   0.15 seconds
            Process time |  65.93 seconds
       Post-process time |   0.00 seconds
              Total time | 192.66 seconds
               User time |   0.00 seconds
       ----------------------------------

      ================================================================================

      Integrating reflections

       Split 4731 reflections overlapping job boundaries

      Processing reflections in the following blocks of images:

       block_size: 34 frames

       --------------------------------------------------------------------------
        # | Group | Frame From | Frame To | Angle From | Angle To | # Reflections
       --------------------------------------------------------------------------
        0 |     0 |          0 |       34 |       82.0 |     87.1 |         17808
        1 |     0 |         17 |       51 |      84.55 |    89.65 |         11824
        2 |     0 |         34 |       68 |       87.1 |     92.2 |         11918
        3 |     0 |         51 |       85 |      89.65 |    94.75 |         11891
        4 |     0 |         68 |      102 |       92.2 |     97.3 |         11806
        5 |     0 |         85 |      119 |      94.75 |    99.85 |         11888
        6 |     0 |        102 |      136 |       97.3 |    102.4 |         11822
        7 |     0 |        119 |      153 |      99.85 |   104.95 |         11823
        8 |     0 |        136 |      170 |      102.4 |    107.5 |         11999
        9 |     0 |        153 |      187 |     104.95 |   110.05 |         11780
       10 |     0 |        170 |      204 |      107.5 |    112.6 |         11867
       11 |     0 |        187 |      221 |     110.05 |   115.15 |         11892
       12 |     0 |        204 |      238 |      112.6 |    117.7 |         11832
       13 |     0 |        221 |      255 |     115.15 |   120.25 |         11862
       14 |     0 |        238 |      272 |      117.7 |    122.8 |         11917
       15 |     0 |        255 |      289 |     120.25 |   125.35 |         11758
       16 |     0 |        272 |      306 |      122.8 |    127.9 |         11798
       17 |     0 |        289 |      323 |     125.35 |   130.45 |         11960
       18 |     0 |        306 |      340 |      127.9 |    133.0 |         11849
       19 |     0 |        323 |      357 |     130.45 |   135.55 |         11794
       20 |     0 |        340 |      374 |      133.0 |    138.1 |         11843
       21 |     0 |        357 |      391 |     135.55 |   140.65 |         11844
       22 |     0 |        374 |      408 |      138.1 |    143.2 |         11883
       23 |     0 |        391 |      425 |     140.65 |   145.75 |         11843
       24 |     0 |        408 |      442 |      143.2 |    148.3 |         11813
       25 |     0 |        425 |      459 |     145.75 |   150.85 |         11810
       26 |     0 |        442 |      476 |      148.3 |    153.4 |         11913
       27 |     0 |        459 |      493 |     150.85 |   155.95 |         11891
       28 |     0 |        476 |      510 |      153.4 |    158.5 |         11857
       29 |     0 |        493 |      527 |     155.95 |   161.05 |         11112
       30 |     0 |        510 |      540 |      158.5 |    163.0 |         15563
       --------------------------------------------------------------------------

       Using multiprocessing with 4 parallel job(s) and 1 thread(s) per job

       Beginning integration job 0

       Frames: 0 -> 34

       Number of reflections
        Partial:     1376
        Full:        16432
        In ice ring: 0
        Integrate:   17808
        Total:       17808

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       0  [2116]: ****************************************
       1  [2842]: ******************************************************
       2  [3244]: **************************************************************
       3  [3321]: ***************************************************************
       4  [3285]: ***************************************************************
       5  [3288]: ***************************************************************
       6  [3316]: ***************************************************************
       7  [3313]: ***************************************************************
       8  [3345]: ****************************************************************
       9  [3385]: *****************************************************************
       10 [3366]: ****************************************************************
       11 [3356]: ****************************************************************
       12 [3323]: ****************************************************************
       13 [3330]: ****************************************************************
       14 [3321]: ***************************************************************
       15 [3394]: *****************************************************************
       16 [3412]: *****************************************************************
       17 [3456]: ******************************************************************
       18 [3495]: *******************************************************************
       19 [3530]: ********************************************************************
       20 [3471]: ******************************************************************
       21 [3398]: *****************************************************************
       22 [3312]: ***************************************************************
       23 [3201]: *************************************************************
       24 [2567]: *************************************************
       25 [1891]: ************************************
       26 [1153]: **********************
       27 [497 ]: *********
       28 [221 ]: ****
       29 [146 ]: **
       30 [111 ]: **
       31 [96  ]: *
       32 [81  ]: *
       33 [69  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0435027 GB

       Integrated   232 (sum) +   231 (prf) /  256 reflections on image 1
       Integrated   559 (sum) +   556 (prf) /  629 reflections on image 2
       Integrated   599 (sum) +   598 (prf) /  681 reflections on image 3
       Integrated   606 (sum) +   601 (prf) /  673 reflections on image 4
       Integrated   583 (sum) +   579 (prf) /  644 reflections on image 5
       Integrated   570 (sum) +   564 (prf) /  660 reflections on image 6
       Integrated   585 (sum) +   577 (prf) /  669 reflections on image 7
       Integrated   573 (sum) +   568 (prf) /  674 reflections on image 8
       Integrated   596 (sum) +   591 (prf) /  679 reflections on image 9
       Integrated   602 (sum) +   597 (prf) /  693 reflections on image 10
       Integrated   625 (sum) +   619 (prf) /  717 reflections on image 11
       Integrated   616 (sum) +   608 (prf) /  694 reflections on image 12
       Integrated   604 (sum) +   600 (prf) /  702 reflections on image 13
       Integrated   550 (sum) +   547 (prf) /  642 reflections on image 14
       Integrated   591 (sum) +   590 (prf) /  670 reflections on image 15
       Integrated   585 (sum) +   580 (prf) /  676 reflections on image 16
       Integrated   578 (sum) +   573 (prf) /  656 reflections on image 17
       Integrated   586 (sum) +   583 (prf) /  676 reflections on image 18
       Integrated   652 (sum) +   644 (prf) /  733 reflections on image 19
       Integrated   653 (sum) +   650 (prf) /  723 reflections on image 20
       Integrated   639 (sum) +   633 (prf) /  724 reflections on image 21
       Integrated   598 (sum) +   592 (prf) /  676 reflections on image 22
       Integrated   621 (sum) +   614 (prf) /  694 reflections on image 23
       Integrated   610 (sum) +   606 (prf) /  676 reflections on image 24
       Integrated   654 (sum) +   650 (prf) /  738 reflections on image 25
       Integrated   582 (sum) +   580 (prf) /  656 reflections on image 26
       Integrated   237 (sum) +   237 (prf) /  276 reflections on image 27
       Integrated    65 (sum) +    64 (prf) /   75 reflections on image 28
       Integrated    29 (sum) +    29 (prf) /   35 reflections on image 29
       Integrated    14 (sum) +    14 (prf) /   15 reflections on image 30
       Integrated    12 (sum) +    12 (prf) /   15 reflections on image 31
       Integrated    11 (sum) +    11 (prf) /   12 reflections on image 32
       Integrated    59 (sum) +    58 (prf) /   69 reflections on image 33
       Beginning integration job 1

       Frames: 17 -> 51

       Number of reflections
        Partial:     162
        Full:        11662
        In ice ring: 0
        Integrate:   11824
        Total:       11824

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       17 [4   ]:
       18 [19  ]:
       19 [30  ]:
       20 [48  ]:
       21 [64  ]: *
       22 [123 ]: **
       23 [244 ]: ****
       24 [873 ]: *****************
       25 [1596]: ********************************
       26 [2318]: **********************************************
       27 [3011]: ************************************************************
       28 [3297]: ******************************************************************
       29 [3337]: ******************************************************************
       30 [3375]: *******************************************************************
       31 [3309]: ******************************************************************
       32 [3313]: ******************************************************************
       33 [3288]: *****************************************************************
       34 [3376]: *******************************************************************
       35 [3339]: ******************************************************************
       36 [3368]: *******************************************************************
       37 [3390]: ********************************************************************
       38 [3317]: ******************************************************************
       39 [3243]: *****************************************************************
       40 [3105]: **************************************************************
       41 [2489]: *************************************************
       42 [1792]: ***********************************
       43 [1136]: **********************
       44 [442 ]: ********
       45 [197 ]: ***
       46 [135 ]: **
       47 [100 ]: **
       48 [80  ]: *
       49 [72  ]: *
       50 [63  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0410384 GB

       Integrated   381 (sum) +   379 (prf) /  439 reflections on image 27
       Integrated   549 (sum) +   546 (prf) /  642 reflections on image 28
       Integrated   581 (sum) +   577 (prf) /  656 reflections on image 29
       Integrated   603 (sum) +   598 (prf) /  693 reflections on image 30
       Integrated   608 (sum) +   605 (prf) /  689 reflections on image 31
       Integrated   626 (sum) +   623 (prf) /  706 reflections on image 32
       Integrated   618 (sum) +   591 (prf) /  711 reflections on image 33
       Integrated   588 (sum) +   578 (prf) /  669 reflections on image 34
       Integrated   587 (sum) +   580 (prf) /  685 reflections on image 35
       Integrated   580 (sum) +   566 (prf) /  662 reflections on image 36
       Integrated   628 (sum) +   620 (prf) /  709 reflections on image 37
       Integrated   598 (sum) +   588 (prf) /  678 reflections on image 38
       Integrated   630 (sum) +   625 (prf) /  716 reflections on image 39
       Integrated   618 (sum) +   614 (prf) /  680 reflections on image 40
       Integrated   621 (sum) +   612 (prf) /  697 reflections on image 41
       Integrated   586 (sum) +   581 (prf) /  656 reflections on image 42
       Integrated   601 (sum) +   599 (prf) /  694 reflections on image 43
       Integrated   217 (sum) +   215 (prf) /  245 reflections on image 44
       Integrated    49 (sum) +    49 (prf) /   62 reflections on image 45
       Integrated    27 (sum) +    27 (prf) /   35 reflections on image 46
       Integrated    19 (sum) +    19 (prf) /   20 reflections on image 47
       Integrated     4 (sum) +     4 (prf) /    8 reflections on image 48
       Integrated     8 (sum) +     8 (prf) /    9 reflections on image 49
       Integrated    50 (sum) +    50 (prf) /   63 reflections on image 50
       Beginning integration job 2

       Frames: 34 -> 68

       Number of reflections
        Partial:     137
        Full:        11781
        In ice ring: 0
        Integrate:   11918
        Total:       11918

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       34 [5   ]:
       35 [14  ]:
       36 [22  ]:
       37 [35  ]:
       38 [60  ]: *
       39 [108 ]: **
       40 [246 ]: ****
       41 [833 ]: ****************
       42 [1561]: ******************************
       43 [2299]: *********************************************
       44 [2977]: **********************************************************
       45 [3292]: ****************************************************************
       46 [3366]: ******************************************************************
       47 [3405]: *******************************************************************
       48 [3437]: *******************************************************************
       49 [3444]: ********************************************************************
       50 [3387]: ******************************************************************
       51 [3434]: *******************************************************************
       52 [3377]: ******************************************************************
       53 [3306]: *****************************************************************
       54 [3342]: *****************************************************************
       55 [3313]: *****************************************************************
       56 [3287]: ****************************************************************
       57 [3159]: **************************************************************
       58 [2591]: ***************************************************
       59 [1878]: *************************************
       60 [1152]: **********************
       61 [471 ]: *********
       62 [207 ]: ****
       63 [138 ]: **
       64 [106 ]: **
       65 [80  ]: *
       66 [71  ]: *
       67 [61  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0417551 GB

       Integrated   347 (sum) +   343 (prf) /  385 reflections on image 44
       Integrated   560 (sum) +   556 (prf) /  628 reflections on image 45
       Integrated   609 (sum) +   607 (prf) /  687 reflections on image 46
       Integrated   613 (sum) +   610 (prf) /  690 reflections on image 47
       Integrated   595 (sum) +   591 (prf) /  695 reflections on image 48
       Integrated   599 (sum) +   597 (prf) /  698 reflections on image 49
       Integrated   651 (sum) +   631 (prf) /  732 reflections on image 50
       Integrated   643 (sum) +   630 (prf) /  746 reflections on image 51
       Integrated   639 (sum) +   634 (prf) /  732 reflections on image 52
       Integrated   554 (sum) +   552 (prf) /  618 reflections on image 53
       Integrated   612 (sum) +   606 (prf) /  694 reflections on image 54
       Integrated   600 (sum) +   592 (prf) /  695 reflections on image 55
       Integrated   603 (sum) +   593 (prf) /  675 reflections on image 56
       Integrated   577 (sum) +   574 (prf) /  652 reflections on image 57
       Integrated   615 (sum) +   611 (prf) /  713 reflections on image 58
       Integrated   651 (sum) +   647 (prf) /  726 reflections on image 59
       Integrated   607 (sum) +   603 (prf) /  681 reflections on image 60
       Integrated   227 (sum) +   222 (prf) /  264 reflections on image 61
       Integrated    59 (sum) +    58 (prf) /   69 reflections on image 62
       Integrated    29 (sum) +    29 (prf) /   32 reflections on image 63
       Integrated    23 (sum) +    23 (prf) /   26 reflections on image 64
       Integrated     6 (sum) +     6 (prf) /    9 reflections on image 65
       Integrated     9 (sum) +     9 (prf) /   10 reflections on image 66
       Integrated    51 (sum) +    51 (prf) /   61 reflections on image 67
       Beginning integration job 3

       Frames: 51 -> 85

       Number of reflections
        Partial:     142
        Full:        11749
        In ice ring: 0
        Integrate:   11891
        Total:       11891

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       51 [5   ]:
       52 [13  ]:
       53 [20  ]:
       54 [45  ]:
       55 [73  ]: *
       56 [118 ]: **
       57 [237 ]: ****
       58 [871 ]: ****************
       59 [1554]: ******************************
       60 [2224]: *******************************************
       61 [2907]: ********************************************************
       62 [3124]: ************************************************************
       63 [3176]: *************************************************************
       64 [3246]: ***************************************************************
       65 [3321]: ****************************************************************
       66 [3399]: ******************************************************************
       67 [3453]: *******************************************************************
       68 [3502]: ********************************************************************
       69 [3468]: *******************************************************************
       70 [3473]: *******************************************************************
       71 [3483]: *******************************************************************
       72 [3450]: ******************************************************************
       73 [3330]: ****************************************************************
       74 [3162]: *************************************************************
       75 [2505]: ************************************************
       76 [1823]: ***********************************
       77 [1169]: **********************
       78 [471 ]: *********
       79 [205 ]: ***
       80 [139 ]: **
       81 [98  ]: *
       82 [74  ]: *
       83 [64  ]: *
       84 [57  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0421338 GB

       Integrated   382 (sum) +   380 (prf) /  440 reflections on image 61
       Integrated   551 (sum) +   549 (prf) /  624 reflections on image 62
       Integrated   555 (sum) +   551 (prf) /  637 reflections on image 63
       Integrated   559 (sum) +   557 (prf) /  638 reflections on image 64
       Integrated   546 (sum) +   545 (prf) /  624 reflections on image 65
       Integrated   590 (sum) +   583 (prf) /  663 reflections on image 66
       Integrated   655 (sum) +   626 (prf) /  742 reflections on image 67
       Integrated   661 (sum) +   649 (prf) /  729 reflections on image 68
       Integrated   619 (sum) +   610 (prf) /  703 reflections on image 69
       Integrated   611 (sum) +   607 (prf) /  709 reflections on image 70
       Integrated   630 (sum) +   625 (prf) /  721 reflections on image 71
       Integrated   623 (sum) +   617 (prf) /  721 reflections on image 72
       Integrated   613 (sum) +   602 (prf) /  696 reflections on image 73
       Integrated   644 (sum) +   639 (prf) /  739 reflections on image 74
       Integrated   596 (sum) +   591 (prf) /  682 reflections on image 75
       Integrated   559 (sum) +   556 (prf) /  654 reflections on image 76
       Integrated   611 (sum) +   607 (prf) /  698 reflections on image 77
       Integrated   226 (sum) +   222 (prf) /  266 reflections on image 78
       Integrated    60 (sum) +    60 (prf) /   66 reflections on image 79
       Integrated    35 (sum) +    35 (prf) /   41 reflections on image 80
       Integrated    21 (sum) +    21 (prf) /   24 reflections on image 81
       Integrated     7 (sum) +     7 (prf) /   10 reflections on image 82
       Integrated     7 (sum) +     7 (prf) /    7 reflections on image 83
       Integrated    47 (sum) +    46 (prf) /   57 reflections on image 84
       Beginning integration job 4

       Frames: 68 -> 102

       Number of reflections
        Partial:     136
        Full:        11670
        In ice ring: 0
        Integrate:   11806
        Total:       11806

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       68  [1   ]:
       69  [18  ]:
       70  [36  ]:
       71  [53  ]:
       72  [80  ]: *
       73  [132 ]: **
       74  [246 ]: ****
       75  [889 ]: ****************
       76  [1518]: ****************************
       77  [2166]: ****************************************
       78  [2861]: *****************************************************
       79  [3137]: **********************************************************
       80  [3220]: ************************************************************
       81  [3281]: *************************************************************
       82  [3299]: *************************************************************
       83  [3359]: ***************************************************************
       84  [3424]: ****************************************************************
       85  [3566]: *******************************************************************
       86  [3551]: ******************************************************************
       87  [3480]: *****************************************************************
       88  [3360]: ***************************************************************
       89  [3315]: **************************************************************
       90  [3303]: **************************************************************
       91  [3172]: ***********************************************************
       92  [2512]: ***********************************************
       93  [1806]: *********************************
       94  [1114]: ********************
       95  [426 ]: ********
       96  [201 ]: ***
       97  [135 ]: **
       98  [96  ]: *
       99  [80  ]: *
       100 [68  ]: *
       101 [58  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0427756 GB

       Integrated   371 (sum) +   369 (prf) /  419 reflections on image 78
       Integrated   519 (sum) +   518 (prf) /  600 reflections on image 79
       Integrated   550 (sum) +   548 (prf) /  621 reflections on image 80
       Integrated   594 (sum) +   588 (prf) /  668 reflections on image 81
       Integrated   560 (sum) +   556 (prf) /  649 reflections on image 82
       Integrated   585 (sum) +   581 (prf) /  680 reflections on image 83
       Integrated   606 (sum) +   584 (prf) /  701 reflections on image 84
       Integrated   628 (sum) +   617 (prf) /  686 reflections on image 85
       Integrated   638 (sum) +   624 (prf) /  720 reflections on image 86
       Integrated   701 (sum) +   697 (prf) /  781 reflections on image 87
       Integrated   614 (sum) +   605 (prf) /  705 reflections on image 88
       Integrated   602 (sum) +   596 (prf) /  683 reflections on image 89
       Integrated   584 (sum) +   581 (prf) /  667 reflections on image 90
       Integrated   632 (sum) +   629 (prf) /  714 reflections on image 91
       Integrated   617 (sum) +   611 (prf) /  706 reflections on image 92
       Integrated   613 (sum) +   606 (prf) /  692 reflections on image 93
       Integrated   617 (sum) +   613 (prf) /  688 reflections on image 94
       Integrated   196 (sum) +   194 (prf) /  225 reflections on image 95
       Integrated    57 (sum) +    57 (prf) /   66 reflections on image 96
       Integrated    34 (sum) +    34 (prf) /   39 reflections on image 97
       Integrated    10 (sum) +    10 (prf) /   16 reflections on image 98
       Integrated    11 (sum) +    11 (prf) /   12 reflections on image 99
       Integrated     8 (sum) +     8 (prf) /   10 reflections on image 100
       Integrated    45 (sum) +    45 (prf) /   58 reflections on image 101
       Beginning integration job 5

       Frames: 85 -> 119

       Number of reflections
        Partial:     140
        Full:        11748
        In ice ring: 0
        Integrate:   11888
        Total:       11888

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       85  [2   ]:
       86  [11  ]:
       87  [24  ]:
       88  [37  ]:
       89  [62  ]: *
       90  [107 ]: **
       91  [251 ]: ****
       92  [854 ]: ****************
       93  [1575]: ******************************
       94  [2266]: *******************************************
       95  [2948]: ********************************************************
       96  [3211]: *************************************************************
       97  [3260]: **************************************************************
       98  [3273]: ***************************************************************
       99  [3334]: ****************************************************************
       100 [3406]: *****************************************************************
       101 [3405]: *****************************************************************
       102 [3474]: *******************************************************************
       103 [3418]: *****************************************************************
       104 [3418]: *****************************************************************
       105 [3455]: ******************************************************************
       106 [3388]: *****************************************************************
       107 [3314]: ***************************************************************
       108 [3170]: *************************************************************
       109 [2521]: ************************************************
       110 [1865]: ***********************************
       111 [1175]: **********************
       112 [465 ]: ********
       113 [228 ]: ****
       114 [140 ]: **
       115 [103 ]: *
       116 [92  ]: *
       117 [76  ]: *
       118 [64  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0423756 GB

       Integrated   370 (sum) +   365 (prf) /  412 reflections on image 95
       Integrated   570 (sum) +   568 (prf) /  658 reflections on image 96
       Integrated   571 (sum) +   569 (prf) /  654 reflections on image 97
       Integrated   599 (sum) +   596 (prf) /  666 reflections on image 98
       Integrated   565 (sum) +   560 (prf) /  649 reflections on image 99
       Integrated   608 (sum) +   603 (prf) /  685 reflections on image 100
       Integrated   624 (sum) +   603 (prf) /  712 reflections on image 101
       Integrated   633 (sum) +   625 (prf) /  722 reflections on image 102
       Integrated   613 (sum) +   603 (prf) /  710 reflections on image 103
       Integrated   603 (sum) +   596 (prf) /  676 reflections on image 104
       Integrated   590 (sum) +   582 (prf) /  681 reflections on image 105
       Integrated   621 (sum) +   615 (prf) /  715 reflections on image 106
       Integrated   608 (sum) +   604 (prf) /  706 reflections on image 107
       Integrated   635 (sum) +   632 (prf) /  721 reflections on image 108
       Integrated   581 (sum) +   576 (prf) /  656 reflections on image 109
       Integrated   591 (sum) +   586 (prf) /  690 reflections on image 110
       Integrated   613 (sum) +   611 (prf) /  710 reflections on image 111
       Integrated   211 (sum) +   209 (prf) /  237 reflections on image 112
       Integrated    75 (sum) +    74 (prf) /   88 reflections on image 113
       Integrated    33 (sum) +    33 (prf) /   37 reflections on image 114
       Integrated     9 (sum) +     9 (prf) /   11 reflections on image 115
       Integrated    16 (sum) +    15 (prf) /   16 reflections on image 116
       Integrated    10 (sum) +    10 (prf) /   12 reflections on image 117
       Integrated    51 (sum) +    51 (prf) /   64 reflections on image 118
       Beginning integration job 6

       Frames: 102 -> 136

       Number of reflections
        Partial:     130
        Full:        11692
        In ice ring: 0
        Integrate:   11822
        Total:       11822

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       102 [7   ]:
       103 [17  ]:
       104 [25  ]:
       105 [41  ]:
       106 [72  ]: *
       107 [108 ]: **
       108 [223 ]: ****
       109 [855 ]: ****************
       110 [1581]: ******************************
       111 [2224]: *******************************************
       112 [2898]: ********************************************************
       113 [3183]: *************************************************************
       114 [3179]: *************************************************************
       115 [3321]: ****************************************************************
       116 [3416]: ******************************************************************
       117 [3400]: *****************************************************************
       118 [3422]: ******************************************************************
       119 [3457]: *******************************************************************
       120 [3404]: *****************************************************************
       121 [3360]: *****************************************************************
       122 [3375]: *****************************************************************
       123 [3365]: *****************************************************************
       124 [3257]: ***************************************************************
       125 [3104]: ************************************************************
       126 [2506]: ************************************************
       127 [1792]: **********************************
       128 [1150]: **********************
       129 [471 ]: *********
       130 [217 ]: ****
       131 [146 ]: **
       132 [100 ]: *
       133 [77  ]: *
       134 [62  ]: *
       135 [57  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0417504 GB

       Integrated   372 (sum) +   368 (prf) /  424 reflections on image 112
       Integrated   588 (sum) +   585 (prf) /  652 reflections on image 113
       Integrated   542 (sum) +   538 (prf) /  610 reflections on image 114
       Integrated   574 (sum) +   571 (prf) /  640 reflections on image 115
       Integrated   646 (sum) +   643 (prf) /  722 reflections on image 116
       Integrated   545 (sum) +   540 (prf) /  627 reflections on image 117
       Integrated   666 (sum) +   644 (prf) /  772 reflections on image 118
       Integrated   643 (sum) +   631 (prf) /  733 reflections on image 119
       Integrated   614 (sum) +   602 (prf) /  706 reflections on image 120
       Integrated   594 (sum) +   588 (prf) /  676 reflections on image 121
       Integrated   603 (sum) +   600 (prf) /  695 reflections on image 122
       Integrated   629 (sum) +   620 (prf) /  700 reflections on image 123
       Integrated   623 (sum) +   615 (prf) /  700 reflections on image 124
       Integrated   588 (sum) +   584 (prf) /  659 reflections on image 125
       Integrated   623 (sum) +   619 (prf) /  714 reflections on image 126
       Integrated   563 (sum) +   556 (prf) /  642 reflections on image 127
       Integrated   590 (sum) +   588 (prf) /  679 reflections on image 128
       Integrated   228 (sum) +   224 (prf) /  254 reflections on image 129
       Integrated    64 (sum) +    64 (prf) /   71 reflections on image 130
       Integrated    41 (sum) +    41 (prf) /   46 reflections on image 131
       Integrated    19 (sum) +    18 (prf) /   23 reflections on image 132
       Integrated    14 (sum) +    14 (prf) /   15 reflections on image 133
       Integrated     4 (sum) +     4 (prf) /    5 reflections on image 134
       Integrated    46 (sum) +    46 (prf) /   57 reflections on image 135
       Beginning integration job 7

       Frames: 119 -> 153

       Number of reflections
        Partial:     149
        Full:        11674
        In ice ring: 0
        Integrate:   11823
        Total:       11823

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       119 [7   ]:
       120 [15  ]:
       121 [31  ]:
       122 [51  ]:
       123 [74  ]: *
       124 [129 ]: **
       125 [272 ]: *****
       126 [876 ]: ****************
       127 [1500]: ****************************
       128 [2184]: *****************************************
       129 [2896]: *******************************************************
       130 [3201]: *************************************************************
       131 [3323]: ***************************************************************
       132 [3443]: ******************************************************************
       133 [3479]: ******************************************************************
       134 [3452]: ******************************************************************
       135 [3452]: ******************************************************************
       136 [3493]: *******************************************************************
       137 [3384]: ****************************************************************
       138 [3355]: ****************************************************************
       139 [3338]: ****************************************************************
       140 [3280]: **************************************************************
       141 [3184]: *************************************************************
       142 [3096]: ***********************************************************
       143 [2490]: ***********************************************
       144 [1804]: **********************************
       145 [1176]: **********************
       146 [505 ]: *********
       147 [230 ]: ****
       148 [155 ]: **
       149 [114 ]: **
       150 [87  ]: *
       151 [76  ]: *
       152 [65  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0421562 GB

       Integrated   394 (sum) +   391 (prf) /  438 reflections on image 129
       Integrated   527 (sum) +   523 (prf) /  588 reflections on image 130
       Integrated   523 (sum) +   515 (prf) /  605 reflections on image 131
       Integrated   620 (sum) +   613 (prf) /  697 reflections on image 132
       Integrated   618 (sum) +   613 (prf) /  714 reflections on image 133
       Integrated   637 (sum) +   633 (prf) /  731 reflections on image 134
       Integrated   638 (sum) +   611 (prf) /  742 reflections on image 135
       Integrated   649 (sum) +   641 (prf) /  745 reflections on image 136
       Integrated   582 (sum) +   573 (prf) /  674 reflections on image 137
       Integrated   615 (sum) +   607 (prf) /  696 reflections on image 138
       Integrated   632 (sum) +   625 (prf) /  707 reflections on image 139
       Integrated   582 (sum) +   576 (prf) /  659 reflections on image 140
       Integrated   587 (sum) +   581 (prf) /  660 reflections on image 141
       Integrated   603 (sum) +   596 (prf) /  677 reflections on image 142
       Integrated   601 (sum) +   599 (prf) /  686 reflections on image 143
       Integrated   550 (sum) +   545 (prf) /  628 reflections on image 144
       Integrated   586 (sum) +   580 (prf) /  671 reflections on image 145
       Integrated   248 (sum) +   247 (prf) /  275 reflections on image 146
       Integrated    64 (sum) +    64 (prf) /   75 reflections on image 147
       Integrated    39 (sum) +    39 (prf) /   41 reflections on image 148
       Integrated    23 (sum) +    23 (prf) /   27 reflections on image 149
       Integrated    10 (sum) +    10 (prf) /   11 reflections on image 150
       Integrated    10 (sum) +    10 (prf) /   11 reflections on image 151
       Integrated    52 (sum) +    52 (prf) /   65 reflections on image 152
       Beginning integration job 8

       Frames: 136 -> 170

       Number of reflections
        Partial:     159
        Full:        11840
        In ice ring: 0
        Integrate:   11999
        Total:       11999

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       136 [8   ]:
       137 [16  ]:
       138 [27  ]:
       139 [36  ]:
       140 [61  ]: *
       141 [110 ]: **
       142 [242 ]: ****
       143 [902 ]: *****************
       144 [1609]: *******************************
       145 [2329]: *********************************************
       146 [3042]: ***********************************************************
       147 [3250]: ***************************************************************
       148 [3316]: ****************************************************************
       149 [3305]: ****************************************************************
       150 [3348]: *****************************************************************
       151 [3343]: *****************************************************************
       152 [3336]: *****************************************************************
       153 [3366]: *****************************************************************
       154 [3341]: *****************************************************************
       155 [3344]: *****************************************************************
       156 [3437]: *******************************************************************
       157 [3436]: ******************************************************************
       158 [3432]: ******************************************************************
       159 [3287]: ****************************************************************
       160 [2596]: **************************************************
       161 [1897]: ************************************
       162 [1168]: **********************
       163 [477 ]: *********
       164 [217 ]: ****
       165 [142 ]: **
       166 [110 ]: **
       167 [93  ]: *
       168 [77  ]: *
       169 [69  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0414766 GB

       Integrated   401 (sum) +   397 (prf) /  461 reflections on image 146
       Integrated   567 (sum) +   561 (prf) /  648 reflections on image 147
       Integrated   608 (sum) +   605 (prf) /  688 reflections on image 148
       Integrated   586 (sum) +   582 (prf) /  670 reflections on image 149
       Integrated   607 (sum) +   602 (prf) /  696 reflections on image 150
       Integrated   573 (sum) +   571 (prf) /  651 reflections on image 151
       Integrated   652 (sum) +   619 (prf) /  732 reflections on image 152
       Integrated   611 (sum) +   598 (prf) /  691 reflections on image 153
       Integrated   626 (sum) +   619 (prf) /  714 reflections on image 154
       Integrated   564 (sum) +   558 (prf) /  633 reflections on image 155
       Integrated   606 (sum) +   601 (prf) /  667 reflections on image 156
       Integrated   595 (sum) +   591 (prf) /  678 reflections on image 157
       Integrated   603 (sum) +   597 (prf) /  709 reflections on image 158
       Integrated   680 (sum) +   670 (prf) /  765 reflections on image 159
       Integrated   613 (sum) +   610 (prf) /  699 reflections on image 160
       Integrated   633 (sum) +   630 (prf) /  729 reflections on image 161
       Integrated   605 (sum) +   597 (prf) /  691 reflections on image 162
       Integrated   231 (sum) +   229 (prf) /  260 reflections on image 163
       Integrated    65 (sum) +    65 (prf) /   75 reflections on image 164
       Integrated    29 (sum) +    29 (prf) /   32 reflections on image 165
       Integrated    17 (sum) +    15 (prf) /   17 reflections on image 166
       Integrated    13 (sum) +    13 (prf) /   16 reflections on image 167
       Integrated     7 (sum) +     7 (prf) /    8 reflections on image 168
       Integrated    61 (sum) +    60 (prf) /   69 reflections on image 169
       Beginning integration job 9

       Frames: 153 -> 187

       Number of reflections
        Partial:     148
        Full:        11632
        In ice ring: 0
        Integrate:   11780
        Total:       11780

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       153 [5   ]:
       154 [18  ]:
       155 [32  ]:
       156 [45  ]:
       157 [70  ]: *
       158 [113 ]: **
       159 [228 ]: ****
       160 [858 ]: ****************
       161 [1551]: *****************************
       162 [2213]: ******************************************
       163 [2903]: *******************************************************
       164 [3206]: *************************************************************
       165 [3240]: **************************************************************
       166 [3249]: **************************************************************
       167 [3346]: ****************************************************************
       168 [3318]: ***************************************************************
       169 [3340]: ****************************************************************
       170 [3479]: ******************************************************************
       171 [3482]: *******************************************************************
       172 [3427]: *****************************************************************
       173 [3400]: *****************************************************************
       174 [3326]: ***************************************************************
       175 [3267]: **************************************************************
       176 [3131]: ************************************************************
       177 [2513]: ************************************************
       178 [1828]: ***********************************
       179 [1131]: *********************
       180 [503 ]: *********
       181 [214 ]: ****
       182 [150 ]: **
       183 [108 ]: **
       184 [88  ]: *
       185 [71  ]: *
       186 [63  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0422319 GB

       Integrated   373 (sum) +   372 (prf) /  410 reflections on image 163
       Integrated   553 (sum) +   545 (prf) /  629 reflections on image 164
       Integrated   579 (sum) +   576 (prf) /  659 reflections on image 165
       Integrated   546 (sum) +   542 (prf) /  615 reflections on image 166
       Integrated   617 (sum) +   611 (prf) /  715 reflections on image 167
       Integrated   605 (sum) +   602 (prf) /  679 reflections on image 168
       Integrated   577 (sum) +   549 (prf) /  676 reflections on image 169
       Integrated   609 (sum) +   600 (prf) /  688 reflections on image 170
       Integrated   641 (sum) +   633 (prf) /  727 reflections on image 171
       Integrated   621 (sum) +   615 (prf) /  706 reflections on image 172
       Integrated   629 (sum) +   622 (prf) /  699 reflections on image 173
       Integrated   620 (sum) +   613 (prf) /  698 reflections on image 174
       Integrated   610 (sum) +   606 (prf) /  688 reflections on image 175
       Integrated   600 (sum) +   598 (prf) /  678 reflections on image 176
       Integrated   597 (sum) +   588 (prf) /  685 reflections on image 177
       Integrated   614 (sum) +   609 (prf) /  697 reflections on image 178
       Integrated   556 (sum) +   551 (prf) /  628 reflections on image 179
       Integrated   250 (sum) +   249 (prf) /  289 reflections on image 180
       Integrated    54 (sum) +    54 (prf) /   64 reflections on image 181
       Integrated    37 (sum) +    37 (prf) /   42 reflections on image 182
       Integrated    17 (sum) +    16 (prf) /   20 reflections on image 183
       Integrated    14 (sum) +    14 (prf) /   17 reflections on image 184
       Integrated     6 (sum) +     6 (prf) /    8 reflections on image 185
       Integrated    56 (sum) +    55 (prf) /   63 reflections on image 186
       Beginning integration job 10

       Frames: 170 -> 204

       Number of reflections
        Partial:     133
        Full:        11734
        In ice ring: 0
        Integrate:   11867
        Total:       11867

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       170 [5   ]:
       171 [9   ]:
       172 [20  ]:
       173 [34  ]:
       174 [53  ]: *
       175 [105 ]: **
       176 [252 ]: ****
       177 [841 ]: ****************
       178 [1513]: *****************************
       179 [2250]: ********************************************
       180 [2975]: **********************************************************
       181 [3237]: ***************************************************************
       182 [3361]: *****************************************************************
       183 [3352]: *****************************************************************
       184 [3366]: *****************************************************************
       185 [3405]: ******************************************************************
       186 [3373]: *****************************************************************
       187 [3426]: *******************************************************************
       188 [3400]: ******************************************************************
       189 [3385]: ******************************************************************
       190 [3310]: ****************************************************************
       191 [3258]: ***************************************************************
       192 [3259]: ***************************************************************
       193 [3212]: **************************************************************
       194 [2585]: **************************************************
       195 [1880]: ************************************
       196 [1206]: ***********************
       197 [503 ]: *********
       198 [210 ]: ****
       199 [136 ]: **
       200 [105 ]: **
       201 [78  ]: *
       202 [64  ]: *
       203 [56  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0411785 GB

       Integrated   365 (sum) +   362 (prf) /  426 reflections on image 180
       Integrated   510 (sum) +   508 (prf) /  581 reflections on image 181
       Integrated   587 (sum) +   582 (prf) /  668 reflections on image 182
       Integrated   609 (sum) +   604 (prf) /  696 reflections on image 183
       Integrated   590 (sum) +   586 (prf) /  668 reflections on image 184
       Integrated   653 (sum) +   647 (prf) /  740 reflections on image 185
       Integrated   628 (sum) +   607 (prf) /  713 reflections on image 186
       Integrated   590 (sum) +   582 (prf) /  683 reflections on image 187
       Integrated   594 (sum) +   587 (prf) /  693 reflections on image 188
       Integrated   646 (sum) +   636 (prf) /  721 reflections on image 189
       Integrated   627 (sum) +   618 (prf) /  709 reflections on image 190
       Integrated   559 (sum) +   552 (prf) /  639 reflections on image 191
       Integrated   592 (sum) +   584 (prf) /  654 reflections on image 192
       Integrated   603 (sum) +   597 (prf) /  691 reflections on image 193
       Integrated   630 (sum) +   628 (prf) /  705 reflections on image 194
       Integrated   595 (sum) +   588 (prf) /  674 reflections on image 195
       Integrated   595 (sum) +   591 (prf) /  703 reflections on image 196
       Integrated   263 (sum) +   260 (prf) /  293 reflections on image 197
       Integrated    67 (sum) +    67 (prf) /   74 reflections on image 198
       Integrated    27 (sum) +    27 (prf) /   31 reflections on image 199
       Integrated    22 (sum) +    22 (prf) /   27 reflections on image 200
       Integrated    12 (sum) +    12 (prf) /   14 reflections on image 201
       Integrated     7 (sum) +     7 (prf) /    8 reflections on image 202
       Integrated    44 (sum) +    44 (prf) /   56 reflections on image 203
       Beginning integration job 11

       Frames: 187 -> 221

       Number of reflections
        Partial:     138
        Full:        11754
        In ice ring: 0
        Integrate:   11892
        Total:       11892

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       187 [3   ]:
       188 [14  ]:
       189 [28  ]:
       190 [39  ]:
       191 [65  ]: *
       192 [113 ]: **
       193 [252 ]: ****
       194 [896 ]: *****************
       195 [1581]: ******************************
       196 [2293]: ********************************************
       197 [2987]: *********************************************************
       198 [3226]: **************************************************************
       199 [3271]: ***************************************************************
       200 [3264]: ***************************************************************
       201 [3277]: ***************************************************************
       202 [3313]: ***************************************************************
       203 [3339]: ****************************************************************
       204 [3377]: *****************************************************************
       205 [3471]: *******************************************************************
       206 [3460]: ******************************************************************
       207 [3457]: ******************************************************************
       208 [3380]: *****************************************************************
       209 [3287]: ***************************************************************
       210 [3193]: *************************************************************
       211 [2546]: *************************************************
       212 [1838]: ***********************************
       213 [1163]: **********************
       214 [488 ]: *********
       215 [208 ]: ****
       216 [139 ]: **
       217 [111 ]: **
       218 [95  ]: *
       219 [79  ]: *
       220 [64  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0419318 GB

       Integrated   384 (sum) +   382 (prf) /  436 reflections on image 197
       Integrated   565 (sum) +   561 (prf) /  643 reflections on image 198
       Integrated   596 (sum) +   595 (prf) /  680 reflections on image 199
       Integrated   577 (sum) +   576 (prf) /  649 reflections on image 200
       Integrated   591 (sum) +   588 (prf) /  671 reflections on image 201
       Integrated   581 (sum) +   580 (prf) /  660 reflections on image 202
       Integrated   640 (sum) +   615 (prf) /  734 reflections on image 203
       Integrated   576 (sum) +   566 (prf) /  653 reflections on image 204
       Integrated   602 (sum) +   591 (prf) /  698 reflections on image 205
       Integrated   614 (sum) +   613 (prf) /  685 reflections on image 206
       Integrated   662 (sum) +   652 (prf) /  746 reflections on image 207
       Integrated   616 (sum) +   610 (prf) /  713 reflections on image 208
       Integrated   590 (sum) +   587 (prf) /  665 reflections on image 209
       Integrated   633 (sum) +   629 (prf) /  713 reflections on image 210
       Integrated   591 (sum) +   586 (prf) /  708 reflections on image 211
       Integrated   581 (sum) +   577 (prf) /  675 reflections on image 212
       Integrated   590 (sum) +   586 (prf) /  675 reflections on image 213
       Integrated   254 (sum) +   253 (prf) /  280 reflections on image 214
       Integrated    61 (sum) +    61 (prf) /   69 reflections on image 215
       Integrated    26 (sum) +    25 (prf) /   28 reflections on image 216
       Integrated    13 (sum) +    13 (prf) /   16 reflections on image 217
       Integrated    12 (sum) +    12 (prf) /   16 reflections on image 218
       Integrated    10 (sum) +    10 (prf) /   15 reflections on image 219
       Integrated    52 (sum) +    51 (prf) /   64 reflections on image 220
       Beginning integration job 12

       Frames: 204 -> 238

       Number of reflections
        Partial:     131
        Full:        11701
        In ice ring: 0
        Integrate:   11832
        Total:       11832

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       204 [4   ]:
       205 [20  ]:
       206 [36  ]:
       207 [50  ]:
       208 [69  ]: *
       209 [115 ]: **
       210 [247 ]: ****
       211 [839 ]: ****************
       212 [1543]: ******************************
       213 [2181]: ******************************************
       214 [2869]: ********************************************************
       215 [3237]: ***************************************************************
       216 [3312]: ****************************************************************
       217 [3356]: *****************************************************************
       218 [3423]: *******************************************************************
       219 [3391]: ******************************************************************
       220 [3381]: ******************************************************************
       221 [3423]: *******************************************************************
       222 [3341]: *****************************************************************
       223 [3292]: ****************************************************************
       224 [3293]: ****************************************************************
       225 [3326]: *****************************************************************
       226 [3285]: ****************************************************************
       227 [3185]: **************************************************************
       228 [2634]: ***************************************************
       229 [1863]: ************************************
       230 [1149]: **********************
       231 [488 ]: *********
       232 [205 ]: ****
       233 [132 ]: **
       234 [90  ]: *
       235 [80  ]: *
       236 [59  ]: *
       237 [55  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.041288 GB

       Integrated   354 (sum) +   351 (prf) /  400 reflections on image 214
       Integrated   546 (sum) +   544 (prf) /  604 reflections on image 215
       Integrated   570 (sum) +   564 (prf) /  651 reflections on image 216
       Integrated   565 (sum) +   560 (prf) /  653 reflections on image 217
       Integrated   653 (sum) +   650 (prf) /  732 reflections on image 218
       Integrated   615 (sum) +   607 (prf) /  701 reflections on image 219
       Integrated   628 (sum) +   608 (prf) /  712 reflections on image 220
       Integrated   661 (sum) +   647 (prf) /  753 reflections on image 221
       Integrated   615 (sum) +   607 (prf) /  673 reflections on image 222
       Integrated   586 (sum) +   580 (prf) /  672 reflections on image 223
       Integrated   566 (sum) +   562 (prf) /  667 reflections on image 224
       Integrated   639 (sum) +   635 (prf) /  712 reflections on image 225
       Integrated   574 (sum) +   570 (prf) /  644 reflections on image 226
       Integrated   561 (sum) +   558 (prf) /  624 reflections on image 227
       Integrated   675 (sum) +   668 (prf) /  771 reflections on image 228
       Integrated   630 (sum) +   626 (prf) /  714 reflections on image 229
       Integrated   578 (sum) +   574 (prf) /  661 reflections on image 230
       Integrated   248 (sum) +   244 (prf) /  283 reflections on image 231
       Integrated    66 (sum) +    66 (prf) /   73 reflections on image 232
       Integrated    32 (sum) +    32 (prf) /   42 reflections on image 233
       Integrated     7 (sum) +     7 (prf) /   10 reflections on image 234
       Integrated    21 (sum) +    21 (prf) /   21 reflections on image 235
       Integrated     1 (sum) +     1 (prf) /    4 reflections on image 236
       Integrated    47 (sum) +    46 (prf) /   55 reflections on image 237
       Beginning integration job 13

       Frames: 221 -> 255

       Number of reflections
        Partial:     140
        Full:        11722
        In ice ring: 0
        Integrate:   11862
        Total:       11862

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       221 [6   ]:
       222 [10  ]:
       223 [22  ]:
       224 [41  ]:
       225 [73  ]: *
       226 [126 ]: **
       227 [239 ]: ****
       228 [896 ]: *****************
       229 [1629]: *******************************
       230 [2309]: ********************************************
       231 [2990]: *********************************************************
       232 [3195]: *************************************************************
       233 [3206]: *************************************************************
       234 [3227]: *************************************************************
       235 [3294]: ***************************************************************
       236 [3368]: ****************************************************************
       237 [3404]: *****************************************************************
       238 [3481]: ******************************************************************
       239 [3503]: *******************************************************************
       240 [3482]: ******************************************************************
       241 [3405]: *****************************************************************
       242 [3352]: ****************************************************************
       243 [3210]: *************************************************************
       244 [3103]: ***********************************************************
       245 [2497]: ***********************************************
       246 [1838]: ***********************************
       247 [1184]: **********************
       248 [482 ]: *********
       249 [226 ]: ****
       250 [144 ]: **
       251 [109 ]: **
       252 [88  ]: *
       253 [71  ]: *
       254 [65  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0423468 GB

       Integrated   396 (sum) +   393 (prf) /  447 reflections on image 231
       Integrated   580 (sum) +   577 (prf) /  654 reflections on image 232
       Integrated   572 (sum) +   571 (prf) /  668 reflections on image 233
       Integrated   568 (sum) +   565 (prf) /  642 reflections on image 234
       Integrated   570 (sum) +   568 (prf) /  660 reflections on image 235
       Integrated   555 (sum) +   550 (prf) /  641 reflections on image 236
       Integrated   604 (sum) +   577 (prf) /  707 reflections on image 237
       Integrated   625 (sum) +   614 (prf) /  716 reflections on image 238
       Integrated   639 (sum) +   628 (prf) /  710 reflections on image 239
       Integrated   647 (sum) +   642 (prf) /  734 reflections on image 240
       Integrated   596 (sum) +   588 (prf) /  669 reflections on image 241
       Integrated   685 (sum) +   676 (prf) /  761 reflections on image 242
       Integrated   614 (sum) +   608 (prf) /  682 reflections on image 243
       Integrated   590 (sum) +   585 (prf) /  674 reflections on image 244
       Integrated   575 (sum) +   569 (prf) /  659 reflections on image 245
       Integrated   581 (sum) +   577 (prf) /  654 reflections on image 246
       Integrated   609 (sum) +   603 (prf) /  702 reflections on image 247
       Integrated   224 (sum) +   221 (prf) /  256 reflections on image 248
       Integrated    76 (sum) +    75 (prf) /   82 reflections on image 249
       Integrated    30 (sum) +    30 (prf) /   35 reflections on image 250
       Integrated    21 (sum) +    21 (prf) /   21 reflections on image 251
       Integrated    16 (sum) +    16 (prf) /   17 reflections on image 252
       Integrated     5 (sum) +     5 (prf) /    6 reflections on image 253
       Integrated    57 (sum) +    57 (prf) /   65 reflections on image 254
       Beginning integration job 14

       Frames: 238 -> 272

       Number of reflections
        Partial:     145
        Full:        11772
        In ice ring: 0
        Integrate:   11917
        Total:       11917

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       238 [6   ]:
       239 [13  ]:
       240 [20  ]:
       241 [37  ]:
       242 [60  ]: *
       243 [108 ]: **
       244 [252 ]: ****
       245 [856 ]: ****************
       246 [1541]: ******************************
       247 [2236]: *******************************************
       248 [2938]: *********************************************************
       249 [3254]: ***************************************************************
       250 [3379]: ******************************************************************
       251 [3430]: *******************************************************************
       252 [3379]: ******************************************************************
       253 [3395]: ******************************************************************
       254 [3333]: *****************************************************************
       255 [3276]: ***************************************************************
       256 [3286]: ****************************************************************
       257 [3284]: ****************************************************************
       258 [3342]: *****************************************************************
       259 [3302]: ****************************************************************
       260 [3314]: ****************************************************************
       261 [3282]: ****************************************************************
       262 [2682]: ****************************************************
       263 [1952]: **************************************
       264 [1231]: ************************
       265 [476 ]: *********
       266 [202 ]: ***
       267 [140 ]: **
       268 [99  ]: *
       269 [79  ]: *
       270 [72  ]: *
       271 [64  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0405809 GB

       Integrated   363 (sum) +   360 (prf) /  417 reflections on image 248
       Integrated   537 (sum) +   536 (prf) /  605 reflections on image 249
       Integrated   595 (sum) +   590 (prf) /  679 reflections on image 250
       Integrated   599 (sum) +   591 (prf) /  706 reflections on image 251
       Integrated   582 (sum) +   577 (prf) /  659 reflections on image 252
       Integrated   653 (sum) +   649 (prf) /  728 reflections on image 253
       Integrated   698 (sum) +   675 (prf) /  787 reflections on image 254
       Integrated   593 (sum) +   585 (prf) /  673 reflections on image 255
       Integrated   574 (sum) +   565 (prf) /  641 reflections on image 256
       Integrated   548 (sum) +   541 (prf) /  629 reflections on image 257
       Integrated   631 (sum) +   624 (prf) /  710 reflections on image 258
       Integrated   587 (sum) +   581 (prf) /  676 reflections on image 259
       Integrated   554 (sum) +   551 (prf) /  647 reflections on image 260
       Integrated   604 (sum) +   602 (prf) /  678 reflections on image 261
       Integrated   629 (sum) +   626 (prf) /  730 reflections on image 262
       Integrated   638 (sum) +   628 (prf) /  721 reflections on image 263
       Integrated   671 (sum) +   668 (prf) /  755 reflections on image 264
       Integrated   231 (sum) +   229 (prf) /  274 reflections on image 265
       Integrated    56 (sum) +    55 (prf) /   62 reflections on image 266
       Integrated    38 (sum) +    38 (prf) /   41 reflections on image 267
       Integrated    18 (sum) +    18 (prf) /   20 reflections on image 268
       Integrated     5 (sum) +     5 (prf) /    7 reflections on image 269
       Integrated     7 (sum) +     7 (prf) /    8 reflections on image 270
       Integrated    56 (sum) +    56 (prf) /   64 reflections on image 271
       Beginning integration job 15

       Frames: 255 -> 289

       Number of reflections
        Partial:     130
        Full:        11628
        In ice ring: 0
        Integrate:   11758
        Total:       11758

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       255 [5   ]:
       256 [11  ]:
       257 [19  ]:
       258 [39  ]:
       259 [71  ]: *
       260 [118 ]: **
       261 [244 ]: ****
       262 [882 ]: ****************
       263 [1582]: ******************************
       264 [2193]: ******************************************
       265 [2870]: *******************************************************
       266 [3104]: ***********************************************************
       267 [3177]: *************************************************************
       268 [3242]: **************************************************************
       269 [3267]: **************************************************************
       270 [3343]: ****************************************************************
       271 [3332]: ****************************************************************
       272 [3456]: ******************************************************************
       273 [3469]: ******************************************************************
       274 [3486]: ******************************************************************
       275 [3487]: *******************************************************************
       276 [3437]: ******************************************************************
       277 [3333]: ****************************************************************
       278 [3110]: ***********************************************************
       279 [2444]: **********************************************
       280 [1764]: *********************************
       281 [1112]: *********************
       282 [469 ]: *********
       283 [197 ]: ***
       284 [126 ]: **
       285 [91  ]: *
       286 [70  ]: *
       287 [57  ]: *
       288 [47  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0420577 GB

       Integrated   383 (sum) +   381 (prf) /  436 reflections on image 265
       Integrated   541 (sum) +   540 (prf) /  606 reflections on image 266
       Integrated   548 (sum) +   540 (prf) /  624 reflections on image 267
       Integrated   575 (sum) +   573 (prf) /  649 reflections on image 268
       Integrated   560 (sum) +   552 (prf) /  630 reflections on image 269
       Integrated   612 (sum) +   604 (prf) /  696 reflections on image 270
       Integrated   604 (sum) +   584 (prf) /  689 reflections on image 271
       Integrated   632 (sum) +   618 (prf) /  713 reflections on image 272
       Integrated   616 (sum) +   610 (prf) /  693 reflections on image 273
       Integrated   616 (sum) +   612 (prf) /  684 reflections on image 274
       Integrated   612 (sum) +   606 (prf) /  695 reflections on image 275
       Integrated   667 (sum) +   660 (prf) /  746 reflections on image 276
       Integrated   627 (sum) +   619 (prf) /  724 reflections on image 277
       Integrated   643 (sum) +   636 (prf) /  729 reflections on image 278
       Integrated   595 (sum) +   588 (prf) /  680 reflections on image 279
       Integrated   554 (sum) +   546 (prf) /  652 reflections on image 280
       Integrated   558 (sum) +   552 (prf) /  643 reflections on image 281
       Integrated   241 (sum) +   238 (prf) /  272 reflections on image 282
       Integrated    66 (sum) +    66 (prf) /   71 reflections on image 283
       Integrated    34 (sum) +    34 (prf) /   35 reflections on image 284
       Integrated    18 (sum) +    18 (prf) /   21 reflections on image 285
       Integrated    13 (sum) +    13 (prf) /   13 reflections on image 286
       Integrated     9 (sum) +     9 (prf) /   10 reflections on image 287
       Integrated    36 (sum) +    36 (prf) /   47 reflections on image 288
       Beginning integration job 16

       Frames: 272 -> 306

       Number of reflections
        Partial:     119
        Full:        11679
        In ice ring: 0
        Integrate:   11798
        Total:       11798

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       272 [10  ]:
       273 [17  ]:
       274 [31  ]:
       275 [49  ]:
       276 [71  ]: *
       277 [112 ]: **
       278 [240 ]: ****
       279 [873 ]: ****************
       280 [1575]: ******************************
       281 [2270]: *******************************************
       282 [2970]: *********************************************************
       283 [3230]: **************************************************************
       284 [3328]: ****************************************************************
       285 [3359]: *****************************************************************
       286 [3407]: *****************************************************************
       287 [3459]: *******************************************************************
       288 [3440]: ******************************************************************
       289 [3446]: ******************************************************************
       290 [3346]: ****************************************************************
       291 [3279]: ***************************************************************
       292 [3276]: ***************************************************************
       293 [3283]: ***************************************************************
       294 [3241]: **************************************************************
       295 [3104]: ************************************************************
       296 [2483]: ************************************************
       297 [1805]: **********************************
       298 [1120]: *********************
       299 [447 ]: ********
       300 [220 ]: ****
       301 [136 ]: **
       302 [100 ]: *
       303 [74  ]: *
       304 [60  ]: *
       305 [51  ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0419309 GB

       Integrated   372 (sum) +   371 (prf) /  420 reflections on image 282
       Integrated   557 (sum) +   553 (prf) /  633 reflections on image 283
       Integrated   586 (sum) +   582 (prf) /  672 reflections on image 284
       Integrated   589 (sum) +   585 (prf) /  675 reflections on image 285
       Integrated   607 (sum) +   604 (prf) /  678 reflections on image 286
       Integrated   600 (sum) +   597 (prf) /  691 reflections on image 287
       Integrated   648 (sum) +   634 (prf) /  731 reflections on image 288
       Integrated   661 (sum) +   655 (prf) /  737 reflections on image 289
       Integrated   652 (sum) +   642 (prf) /  729 reflections on image 290
       Integrated   576 (sum) +   570 (prf) /  667 reflections on image 291
       Integrated   591 (sum) +   581 (prf) /  661 reflections on image 292
       Integrated   593 (sum) +   589 (prf) /  673 reflections on image 293
       Integrated   596 (sum) +   589 (prf) /  672 reflections on image 294
       Integrated   590 (sum) +   586 (prf) /  676 reflections on image 295
       Integrated   581 (sum) +   577 (prf) /  678 reflections on image 296
       Integrated   587 (sum) +   583 (prf) /  685 reflections on image 297
       Integrated   574 (sum) +   571 (prf) /  673 reflections on image 298
       Integrated   199 (sum) +   196 (prf) /  227 reflections on image 299
       Integrated    76 (sum) +    76 (prf) /   84 reflections on image 300
       Integrated    33 (sum) +    33 (prf) /   36 reflections on image 301
       Integrated    24 (sum) +    24 (prf) /   26 reflections on image 302
       Integrated    12 (sum) +    12 (prf) /   14 reflections on image 303
       Integrated     6 (sum) +     6 (prf) /    9 reflections on image 304
       Integrated    40 (sum) +    40 (prf) /   51 reflections on image 305
       Beginning integration job 17

       Frames: 289 -> 323

       Number of reflections
        Partial:     132
        Full:        11828
        In ice ring: 0
        Integrate:   11960
        Total:       11960

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       289 [7   ]:
       290 [18  ]:
       291 [32  ]:
       292 [54  ]: *
       293 [77  ]: *
       294 [126 ]: **
       295 [259 ]: *****
       296 [853 ]: ****************
       297 [1567]: ******************************
       298 [2276]: ********************************************
       299 [3027]: ***********************************************************
       300 [3369]: *****************************************************************
       301 [3434]: *******************************************************************
       302 [3417]: ******************************************************************
       303 [3343]: *****************************************************************
       304 [3324]: ****************************************************************
       305 [3265]: ***************************************************************
       306 [3344]: *****************************************************************
       307 [3351]: *****************************************************************
       308 [3328]: ****************************************************************
       309 [3284]: ****************************************************************
       310 [3285]: ****************************************************************
       311 [3323]: ****************************************************************
       312 [3235]: ***************************************************************
       313 [2648]: ***************************************************
       314 [1920]: *************************************
       315 [1229]: ***********************
       316 [513 ]: **********
       317 [225 ]: ****
       318 [134 ]: **
       319 [100 ]: *
       320 [78  ]: *
       321 [63  ]: *
       322 [53  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0408243 GB

       Integrated   360 (sum) +   358 (prf) /  399 reflections on image 299
       Integrated   573 (sum) +   567 (prf) /  663 reflections on image 300
       Integrated   610 (sum) +   601 (prf) /  693 reflections on image 301
       Integrated   632 (sum) +   625 (prf) /  730 reflections on image 302
       Integrated   623 (sum) +   621 (prf) /  699 reflections on image 303
       Integrated   652 (sum) +   646 (prf) /  734 reflections on image 304
       Integrated   602 (sum) +   579 (prf) /  681 reflections on image 305
       Integrated   597 (sum) +   589 (prf) /  674 reflections on image 306
       Integrated   591 (sum) +   579 (prf) /  672 reflections on image 307
       Integrated   610 (sum) +   609 (prf) /  694 reflections on image 308
       Integrated   594 (sum) +   589 (prf) /  683 reflections on image 309
       Integrated   560 (sum) +   557 (prf) /  645 reflections on image 310
       Integrated   614 (sum) +   609 (prf) /  693 reflections on image 311
       Integrated   574 (sum) +   566 (prf) /  652 reflections on image 312
       Integrated   662 (sum) +   656 (prf) /  728 reflections on image 313
       Integrated   610 (sum) +   608 (prf) /  691 reflections on image 314
       Integrated   622 (sum) +   620 (prf) /  716 reflections on image 315
       Integrated   258 (sum) +   256 (prf) /  288 reflections on image 316
       Integrated    88 (sum) +    87 (prf) /   91 reflections on image 317
       Integrated    30 (sum) +    30 (prf) /   34 reflections on image 318
       Integrated    18 (sum) +    18 (prf) /   22 reflections on image 319
       Integrated    13 (sum) +    13 (prf) /   15 reflections on image 320
       Integrated     7 (sum) +     7 (prf) /   10 reflections on image 321
       Integrated    46 (sum) +    46 (prf) /   53 reflections on image 322
       Beginning integration job 18

       Frames: 306 -> 340

       Number of reflections
        Partial:     128
        Full:        11721
        In ice ring: 0
        Integrate:   11849
        Total:       11849

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       306 [2   ]:
       307 [14  ]:
       308 [25  ]:
       309 [42  ]:
       310 [62  ]: *
       311 [110 ]: **
       312 [230 ]: ****
       313 [857 ]: ****************
       314 [1536]: *****************************
       315 [2260]: ********************************************
       316 [2954]: *********************************************************
       317 [3223]: **************************************************************
       318 [3216]: **************************************************************
       319 [3275]: ***************************************************************
       320 [3254]: ***************************************************************
       321 [3292]: ****************************************************************
       322 [3302]: ****************************************************************
       323 [3402]: ******************************************************************
       324 [3433]: *******************************************************************
       325 [3399]: ******************************************************************
       326 [3360]: *****************************************************************
       327 [3393]: ******************************************************************
       328 [3354]: *****************************************************************
       329 [3273]: ***************************************************************
       330 [2615]: ***************************************************
       331 [1904]: *************************************
       332 [1170]: **********************
       333 [494 ]: *********
       334 [206 ]: ****
       335 [126 ]: **
       336 [90  ]: *
       337 [72  ]: *
       338 [62  ]: *
       339 [52  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0413676 GB

       Integrated   367 (sum) +   365 (prf) /  415 reflections on image 316
       Integrated   561 (sum) +   556 (prf) /  637 reflections on image 317
       Integrated   557 (sum) +   555 (prf) /  648 reflections on image 318
       Integrated   601 (sum) +   595 (prf) /  681 reflections on image 319
       Integrated   586 (sum) +   583 (prf) /  662 reflections on image 320
       Integrated   593 (sum) +   587 (prf) /  666 reflections on image 321
       Integrated   613 (sum) +   589 (prf) /  701 reflections on image 322
       Integrated   574 (sum) +   567 (prf) /  661 reflections on image 323
       Integrated   618 (sum) +   610 (prf) /  712 reflections on image 324
       Integrated   617 (sum) +   614 (prf) /  698 reflections on image 325
       Integrated   583 (sum) +   573 (prf) /  677 reflections on image 326
       Integrated   624 (sum) +   619 (prf) /  708 reflections on image 327
       Integrated   595 (sum) +   584 (prf) /  659 reflections on image 328
       Integrated   619 (sum) +   616 (prf) /  709 reflections on image 329
       Integrated   628 (sum) +   626 (prf) /  711 reflections on image 330
       Integrated   648 (sum) +   644 (prf) /  734 reflections on image 331
       Integrated   608 (sum) +   600 (prf) /  676 reflections on image 332
       Integrated   255 (sum) +   252 (prf) /  288 reflections on image 333
       Integrated    70 (sum) +    69 (prf) /   80 reflections on image 334
       Integrated    33 (sum) +    33 (prf) /   36 reflections on image 335
       Integrated    12 (sum) +    12 (prf) /   18 reflections on image 336
       Integrated     9 (sum) +     9 (prf) /   10 reflections on image 337
       Integrated     8 (sum) +     8 (prf) /   10 reflections on image 338
       Integrated    39 (sum) +    39 (prf) /   52 reflections on image 339
       Beginning integration job 19

       Frames: 323 -> 357

       Number of reflections
        Partial:     127
        Full:        11667
        In ice ring: 0
        Integrate:   11794
        Total:       11794

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       323 [7   ]:
       324 [16  ]:
       325 [31  ]:
       326 [41  ]:
       327 [63  ]: *
       328 [109 ]: **
       329 [230 ]: ****
       330 [873 ]: *****************
       331 [1570]: ******************************
       332 [2230]: *******************************************
       333 [2943]: *********************************************************
       334 [3142]: *************************************************************
       335 [3201]: **************************************************************
       336 [3290]: ****************************************************************
       337 [3305]: ****************************************************************
       338 [3342]: *****************************************************************
       339 [3394]: ******************************************************************
       340 [3406]: ******************************************************************
       341 [3413]: *******************************************************************
       342 [3368]: ******************************************************************
       343 [3339]: *****************************************************************
       344 [3383]: ******************************************************************
       345 [3306]: ****************************************************************
       346 [3152]: *************************************************************
       347 [2522]: *************************************************
       348 [1841]: ************************************
       349 [1181]: ***********************
       350 [513 ]: **********
       351 [207 ]: ****
       352 [144 ]: **
       353 [102 ]: **
       354 [83  ]: *
       355 [68  ]: *
       356 [57  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0415398 GB

       Integrated   394 (sum) +   391 (prf) /  443 reflections on image 333
       Integrated   545 (sum) +   543 (prf) /  619 reflections on image 334
       Integrated   549 (sum) +   548 (prf) /  625 reflections on image 335
       Integrated   582 (sum) +   579 (prf) /  658 reflections on image 336
       Integrated   584 (sum) +   577 (prf) /  684 reflections on image 337
       Integrated   565 (sum) +   562 (prf) /  633 reflections on image 338
       Integrated   641 (sum) +   615 (prf) /  739 reflections on image 339
       Integrated   606 (sum) +   600 (prf) /  690 reflections on image 340
       Integrated   610 (sum) +   598 (prf) /  716 reflections on image 341
       Integrated   638 (sum) +   622 (prf) /  723 reflections on image 342
       Integrated   553 (sum) +   545 (prf) /  635 reflections on image 343
       Integrated   599 (sum) +   592 (prf) /  692 reflections on image 344
       Integrated   636 (sum) +   631 (prf) /  718 reflections on image 345
       Integrated   608 (sum) +   602 (prf) /  697 reflections on image 346
       Integrated   614 (sum) +   612 (prf) /  681 reflections on image 347
       Integrated   584 (sum) +   582 (prf) /  660 reflections on image 348
       Integrated   585 (sum) +   581 (prf) /  668 reflections on image 349
       Integrated   263 (sum) +   262 (prf) /  306 reflections on image 350
       Integrated    54 (sum) +    53 (prf) /   63 reflections on image 351
       Integrated    36 (sum) +    36 (prf) /   42 reflections on image 352
       Integrated    18 (sum) +    18 (prf) /   19 reflections on image 353
       Integrated    12 (sum) +    12 (prf) /   15 reflections on image 354
       Integrated     8 (sum) +     8 (prf) /   11 reflections on image 355
       Integrated    45 (sum) +    44 (prf) /   57 reflections on image 356
       Beginning integration job 20

       Frames: 340 -> 374

       Number of reflections
        Partial:     121
        Full:        11722
        In ice ring: 0
        Integrate:   11843
        Total:       11843

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       340 [2   ]:
       341 [10  ]:
       342 [27  ]:
       343 [47  ]:
       344 [67  ]: *
       345 [128 ]: **
       346 [253 ]: *****
       347 [907 ]: ******************
       348 [1653]: *********************************
       349 [2337]: **********************************************
       350 [3022]: ************************************************************
       351 [3238]: ****************************************************************
       352 [3203]: ****************************************************************
       353 [3292]: *****************************************************************
       354 [3238]: ****************************************************************
       355 [3296]: ******************************************************************
       356 [3291]: *****************************************************************
       357 [3344]: *******************************************************************
       358 [3337]: ******************************************************************
       359 [3330]: ******************************************************************
       360 [3312]: ******************************************************************
       361 [3301]: ******************************************************************
       362 [3318]: ******************************************************************
       363 [3228]: ****************************************************************
       364 [2560]: ***************************************************
       365 [1887]: *************************************
       366 [1163]: ***********************
       367 [477 ]: *********
       368 [205 ]: ****
       369 [131 ]: **
       370 [96  ]: *
       371 [76  ]: *
       372 [65  ]: *
       373 [50  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0409209 GB

       Integrated   402 (sum) +   399 (prf) /  456 reflections on image 350
       Integrated   604 (sum) +   600 (prf) /  689 reflections on image 351
       Integrated   558 (sum) +   556 (prf) /  624 reflections on image 352
       Integrated   595 (sum) +   588 (prf) /  684 reflections on image 353
       Integrated   588 (sum) +   581 (prf) /  666 reflections on image 354
       Integrated   584 (sum) +   580 (prf) /  677 reflections on image 355
       Integrated   606 (sum) +   583 (prf) /  701 reflections on image 356
       Integrated   605 (sum) +   599 (prf) /  671 reflections on image 357
       Integrated   593 (sum) +   585 (prf) /  670 reflections on image 358
       Integrated   614 (sum) +   608 (prf) /  680 reflections on image 359
       Integrated   596 (sum) +   587 (prf) /  688 reflections on image 360
       Integrated   594 (sum) +   586 (prf) /  662 reflections on image 361
       Integrated   603 (sum) +   599 (prf) /  675 reflections on image 362
       Integrated   644 (sum) +   639 (prf) /  740 reflections on image 363
       Integrated   592 (sum) +   587 (prf) /  673 reflections on image 364
       Integrated   640 (sum) +   632 (prf) /  724 reflections on image 365
       Integrated   604 (sum) +   598 (prf) /  686 reflections on image 366
       Integrated   243 (sum) +   241 (prf) /  272 reflections on image 367
       Integrated    68 (sum) +    68 (prf) /   74 reflections on image 368
       Integrated    31 (sum) +    31 (prf) /   35 reflections on image 369
       Integrated    17 (sum) +    17 (prf) /   20 reflections on image 370
       Integrated     9 (sum) +     9 (prf) /   11 reflections on image 371
       Integrated    11 (sum) +    11 (prf) /   15 reflections on image 372
       Integrated    41 (sum) +    41 (prf) /   50 reflections on image 373
       Beginning integration job 21

       Frames: 357 -> 391

       Number of reflections
        Partial:     132
        Full:        11712
        In ice ring: 0
        Integrate:   11844
        Total:       11844

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       357 [7   ]:
       358 [16  ]:
       359 [26  ]:
       360 [50  ]:
       361 [78  ]: *
       362 [132 ]: **
       363 [254 ]: *****
       364 [901 ]: *****************
       365 [1563]: ******************************
       366 [2310]: *********************************************
       367 [3025]: ***********************************************************
       368 [3277]: ****************************************************************
       369 [3311]: *****************************************************************
       370 [3293]: *****************************************************************
       371 [3262]: ****************************************************************
       372 [3283]: ****************************************************************
       373 [3269]: ****************************************************************
       374 [3365]: ******************************************************************
       375 [3388]: ******************************************************************
       376 [3389]: *******************************************************************
       377 [3378]: ******************************************************************
       378 [3337]: *****************************************************************
       379 [3303]: *****************************************************************
       380 [3170]: **************************************************************
       381 [2541]: **************************************************
       382 [1861]: ************************************
       383 [1172]: ***********************
       384 [499 ]: *********
       385 [197 ]: ***
       386 [125 ]: **
       387 [98  ]: *
       388 [84  ]: *
       389 [68  ]: *
       390 [54  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.041025 GB

       Integrated   384 (sum) +   381 (prf) /  438 reflections on image 367
       Integrated   542 (sum) +   534 (prf) /  617 reflections on image 368
       Integrated   620 (sum) +   617 (prf) /  703 reflections on image 369
       Integrated   612 (sum) +   609 (prf) /  689 reflections on image 370
       Integrated   595 (sum) +   591 (prf) /  670 reflections on image 371
       Integrated   619 (sum) +   611 (prf) /  712 reflections on image 372
       Integrated   599 (sum) +   573 (prf) /  677 reflections on image 373
       Integrated   572 (sum) +   563 (prf) /  648 reflections on image 374
       Integrated   612 (sum) +   605 (prf) /  694 reflections on image 375
       Integrated   610 (sum) +   607 (prf) /  700 reflections on image 376
       Integrated   579 (sum) +   573 (prf) /  668 reflections on image 377
       Integrated   612 (sum) +   606 (prf) /  694 reflections on image 378
       Integrated   616 (sum) +   611 (prf) /  699 reflections on image 379
       Integrated   601 (sum) +   598 (prf) /  694 reflections on image 380
       Integrated   595 (sum) +   589 (prf) /  680 reflections on image 381
       Integrated   616 (sum) +   613 (prf) /  689 reflections on image 382
       Integrated   597 (sum) +   596 (prf) /  673 reflections on image 383
       Integrated   255 (sum) +   254 (prf) /  302 reflections on image 384
       Integrated    64 (sum) +    64 (prf) /   72 reflections on image 385
       Integrated    23 (sum) +    23 (prf) /   27 reflections on image 386
       Integrated    14 (sum) +    14 (prf) /   14 reflections on image 387
       Integrated    14 (sum) +    13 (prf) /   16 reflections on image 388
       Integrated    11 (sum) +    11 (prf) /   14 reflections on image 389
       Integrated    44 (sum) +    43 (prf) /   54 reflections on image 390
       Beginning integration job 22

       Frames: 374 -> 408

       Number of reflections
        Partial:     141
        Full:        11742
        In ice ring: 0
        Integrate:   11883
        Total:       11883

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       374 [6   ]:
       375 [14  ]:
       376 [25  ]:
       377 [36  ]:
       378 [58  ]: *
       379 [113 ]: **
       380 [239 ]: ****
       381 [817 ]: ***************
       382 [1501]: ****************************
       383 [2221]: ******************************************
       384 [2947]: ********************************************************
       385 [3225]: **************************************************************
       386 [3358]: ****************************************************************
       387 [3402]: *****************************************************************
       388 [3413]: *****************************************************************
       389 [3443]: ******************************************************************
       390 [3439]: ******************************************************************
       391 [3470]: *******************************************************************
       392 [3369]: *****************************************************************
       393 [3359]: ****************************************************************
       394 [3321]: ****************************************************************
       395 [3284]: ***************************************************************
       396 [3257]: **************************************************************
       397 [3150]: ************************************************************
       398 [2553]: *************************************************
       399 [1866]: ************************************
       400 [1165]: **********************
       401 [475 ]: *********
       402 [209 ]: ****
       403 [149 ]: **
       404 [108 ]: **
       405 [87  ]: *
       406 [70  ]: *
       407 [63  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0415927 GB

       Integrated   345 (sum) +   344 (prf) /  393 reflections on image 384
       Integrated   527 (sum) +   523 (prf) /  601 reflections on image 385
       Integrated   596 (sum) +   587 (prf) /  675 reflections on image 386
       Integrated   622 (sum) +   620 (prf) /  703 reflections on image 387
       Integrated   589 (sum) +   581 (prf) /  670 reflections on image 388
       Integrated   636 (sum) +   633 (prf) /  699 reflections on image 389
       Integrated   654 (sum) +   633 (prf) /  754 reflections on image 390
       Integrated   662 (sum) +   652 (prf) /  739 reflections on image 391
       Integrated   591 (sum) +   589 (prf) /  676 reflections on image 392
       Integrated   606 (sum) +   602 (prf) /  700 reflections on image 393
       Integrated   603 (sum) +   599 (prf) /  683 reflections on image 394
       Integrated   585 (sum) +   581 (prf) /  682 reflections on image 395
       Integrated   602 (sum) +   596 (prf) /  679 reflections on image 396
       Integrated   585 (sum) +   577 (prf) /  676 reflections on image 397
       Integrated   592 (sum) +   590 (prf) /  687 reflections on image 398
       Integrated   615 (sum) +   605 (prf) /  701 reflections on image 399
       Integrated   607 (sum) +   600 (prf) /  690 reflections on image 400
       Integrated   238 (sum) +   236 (prf) /  266 reflections on image 401
       Integrated    55 (sum) +    54 (prf) /   60 reflections on image 402
       Integrated    36 (sum) +    36 (prf) /   41 reflections on image 403
       Integrated    21 (sum) +    20 (prf) /   21 reflections on image 404
       Integrated    15 (sum) +    15 (prf) /   17 reflections on image 405
       Integrated     7 (sum) +     7 (prf) /    7 reflections on image 406
       Integrated    52 (sum) +    52 (prf) /   63 reflections on image 407
       Beginning integration job 23

       Frames: 391 -> 425

       Number of reflections
        Partial:     129
        Full:        11714
        In ice ring: 0
        Integrate:   11843
        Total:       11843

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       392 [11  ]:
       393 [25  ]:
       394 [35  ]:
       395 [63  ]: *
       396 [116 ]: **
       397 [247 ]: ****
       398 [854 ]: ****************
       399 [1552]: ******************************
       400 [2214]: *******************************************
       401 [2852]: *******************************************************
       402 [3138]: *************************************************************
       403 [3238]: ***************************************************************
       404 [3318]: ****************************************************************
       405 [3357]: *****************************************************************
       406 [3362]: *****************************************************************
       407 [3422]: ******************************************************************
       408 [3427]: ******************************************************************
       409 [3384]: ******************************************************************
       410 [3428]: *******************************************************************
       411 [3399]: ******************************************************************
       412 [3353]: *****************************************************************
       413 [3305]: ****************************************************************
       414 [3183]: **************************************************************
       415 [2603]: **************************************************
       416 [1910]: *************************************
       417 [1147]: **********************
       418 [482 ]: *********
       419 [214 ]: ****
       420 [136 ]: **
       421 [103 ]: **
       422 [76  ]: *
       423 [63  ]: *
       424 [53  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0411455 GB

       Integrated   385 (sum) +   384 (prf) /  420 reflections on image 401
       Integrated   531 (sum) +   527 (prf) /  604 reflections on image 402
       Integrated   555 (sum) +   546 (prf) /  636 reflections on image 403
       Integrated   554 (sum) +   553 (prf) /  632 reflections on image 404
       Integrated   593 (sum) +   589 (prf) /  672 reflections on image 405
       Integrated   588 (sum) +   585 (prf) /  675 reflections on image 406
       Integrated   702 (sum) +   688 (prf) /  787 reflections on image 407
       Integrated   618 (sum) +   609 (prf) /  703 reflections on image 408
       Integrated   565 (sum) +   560 (prf) /  649 reflections on image 409
       Integrated   618 (sum) +   614 (prf) /  703 reflections on image 410
       Integrated   626 (sum) +   618 (prf) /  707 reflections on image 411
       Integrated   622 (sum) +   616 (prf) /  707 reflections on image 412
       Integrated   602 (sum) +   598 (prf) /  688 reflections on image 413
       Integrated   568 (sum) +   558 (prf) /  657 reflections on image 414
       Integrated   618 (sum) +   612 (prf) /  693 reflections on image 415
       Integrated   675 (sum) +   668 (prf) /  763 reflections on image 416
       Integrated   560 (sum) +   555 (prf) /  665 reflections on image 417
       Integrated   241 (sum) +   240 (prf) /  268 reflections on image 418
       Integrated    70 (sum) +    70 (prf) /   78 reflections on image 419
       Integrated    29 (sum) +    28 (prf) /   33 reflections on image 420
       Integrated    26 (sum) +    26 (prf) /   27 reflections on image 421
       Integrated    12 (sum) +    12 (prf) /   13 reflections on image 422
       Integrated     9 (sum) +     9 (prf) /   10 reflections on image 423
       Integrated    43 (sum) +    43 (prf) /   53 reflections on image 424
       Beginning integration job 24

       Frames: 408 -> 442

       Number of reflections
        Partial:     134
        Full:        11679
        In ice ring: 0
        Integrate:   11813
        Total:       11813

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       408 [4   ]:
       409 [8   ]:
       410 [17  ]:
       411 [30  ]:
       412 [60  ]: *
       413 [103 ]: **
       414 [247 ]: ****
       415 [875 ]: *****************
       416 [1581]: ******************************
       417 [2296]: ********************************************
       418 [3029]: ***********************************************************
       419 [3291]: ****************************************************************
       420 [3344]: *****************************************************************
       421 [3316]: ****************************************************************
       422 [3275]: ****************************************************************
       423 [3315]: ****************************************************************
       424 [3368]: *****************************************************************
       425 [3423]: *******************************************************************
       426 [3341]: *****************************************************************
       427 [3307]: ****************************************************************
       428 [3311]: ****************************************************************
       429 [3278]: ****************************************************************
       430 [3204]: **************************************************************
       431 [3139]: *************************************************************
       432 [2490]: ************************************************
       433 [1820]: ***********************************
       434 [1192]: ***********************
       435 [470 ]: *********
       436 [224 ]: ****
       437 [143 ]: **
       438 [104 ]: **
       439 [86  ]: *
       440 [68  ]: *
       441 [55  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0411977 GB

       Integrated   369 (sum) +   368 (prf) /  422 reflections on image 418
       Integrated   560 (sum) +   555 (prf) /  635 reflections on image 419
       Integrated   581 (sum) +   576 (prf) /  660 reflections on image 420
       Integrated   642 (sum) +   637 (prf) /  726 reflections on image 421
       Integrated   582 (sum) +   579 (prf) /  666 reflections on image 422
       Integrated   583 (sum) +   580 (prf) /  662 reflections on image 423
       Integrated   645 (sum) +   622 (prf) /  708 reflections on image 424
       Integrated   618 (sum) +   613 (prf) /  726 reflections on image 425
       Integrated   611 (sum) +   601 (prf) /  688 reflections on image 426
       Integrated   598 (sum) +   591 (prf) /  704 reflections on image 427
       Integrated   600 (sum) +   589 (prf) /  671 reflections on image 428
       Integrated   600 (sum) +   592 (prf) /  676 reflections on image 429
       Integrated   577 (sum) +   574 (prf) /  665 reflections on image 430
       Integrated   625 (sum) +   621 (prf) /  714 reflections on image 431
       Integrated   598 (sum) +   594 (prf) /  670 reflections on image 432
       Integrated   557 (sum) +   551 (prf) /  628 reflections on image 433
       Integrated   643 (sum) +   637 (prf) /  722 reflections on image 434
       Integrated   219 (sum) +   216 (prf) /  246 reflections on image 435
       Integrated    72 (sum) +    72 (prf) /   81 reflections on image 436
       Integrated    31 (sum) +    31 (prf) /   39 reflections on image 437
       Integrated    14 (sum) +    14 (prf) /   18 reflections on image 438
       Integrated    15 (sum) +    15 (prf) /   18 reflections on image 439
       Integrated     9 (sum) +     9 (prf) /   13 reflections on image 440
       Integrated    50 (sum) +    50 (prf) /   55 reflections on image 441
       Beginning integration job 25

       Frames: 425 -> 459

       Number of reflections
        Partial:     130
        Full:        11680
        In ice ring: 0
        Integrate:   11810
        Total:       11810

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       425 [10  ]:
       426 [22  ]:
       427 [34  ]:
       428 [46  ]:
       429 [62  ]: *
       430 [107 ]: **
       431 [222 ]: ****
       432 [860 ]: ****************
       433 [1576]: ******************************
       434 [2261]: *******************************************
       435 [2968]: *********************************************************
       436 [3184]: *************************************************************
       437 [3170]: *************************************************************
       438 [3194]: *************************************************************
       439 [3282]: ***************************************************************
       440 [3354]: ****************************************************************
       441 [3404]: *****************************************************************
       442 [3475]: *******************************************************************
       443 [3455]: ******************************************************************
       444 [3377]: *****************************************************************
       445 [3374]: *****************************************************************
       446 [3316]: ***************************************************************
       447 [3212]: *************************************************************
       448 [3138]: ************************************************************
       449 [2534]: ************************************************
       450 [1881]: ************************************
       451 [1156]: **********************
       452 [474 ]: *********
       453 [219 ]: ****
       454 [139 ]: **
       455 [102 ]: *
       456 [83  ]: *
       457 [63  ]: *
       458 [53  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0417944 GB

       Integrated   380 (sum) +   379 (prf) /  426 reflections on image 435
       Integrated   593 (sum) +   591 (prf) /  670 reflections on image 436
       Integrated   573 (sum) +   570 (prf) /  657 reflections on image 437
       Integrated   547 (sum) +   541 (prf) /  627 reflections on image 438
       Integrated   567 (sum) +   564 (prf) /  646 reflections on image 439
       Integrated   600 (sum) +   595 (prf) /  665 reflections on image 440
       Integrated   613 (sum) +   594 (prf) /  702 reflections on image 441
       Integrated   599 (sum) +   584 (prf) /  695 reflections on image 442
       Integrated   663 (sum) +   658 (prf) /  743 reflections on image 443
       Integrated   605 (sum) +   595 (prf) /  711 reflections on image 444
       Integrated   603 (sum) +   594 (prf) /  696 reflections on image 445
       Integrated   595 (sum) +   589 (prf) /  676 reflections on image 446
       Integrated   596 (sum) +   592 (prf) /  681 reflections on image 447
       Integrated   591 (sum) +   586 (prf) /  681 reflections on image 448
       Integrated   576 (sum) +   569 (prf) /  653 reflections on image 449
       Integrated   657 (sum) +   650 (prf) /  725 reflections on image 450
       Integrated   600 (sum) +   595 (prf) /  682 reflections on image 451
       Integrated   222 (sum) +   222 (prf) /  255 reflections on image 452
       Integrated    62 (sum) +    62 (prf) /   80 reflections on image 453
       Integrated    34 (sum) +    34 (prf) /   37 reflections on image 454
       Integrated    17 (sum) +    17 (prf) /   19 reflections on image 455
       Integrated    19 (sum) +    18 (prf) /   20 reflections on image 456
       Integrated    10 (sum) +    10 (prf) /   10 reflections on image 457
       Integrated    43 (sum) +    43 (prf) /   53 reflections on image 458
       Beginning integration job 26

       Frames: 442 -> 476

       Number of reflections
        Partial:     128
        Full:        11785
        In ice ring: 0
        Integrate:   11913
        Total:       11913

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       442 [8   ]:
       443 [17  ]:
       444 [25  ]:
       445 [36  ]:
       446 [62  ]: *
       447 [106 ]: **
       448 [240 ]: ****
       449 [890 ]: *****************
       450 [1610]: *******************************
       451 [2342]: *********************************************
       452 [2995]: **********************************************************
       453 [3213]: **************************************************************
       454 [3234]: ***************************************************************
       455 [3237]: ***************************************************************
       456 [3339]: *****************************************************************
       457 [3318]: ****************************************************************
       458 [3378]: *****************************************************************
       459 [3429]: ******************************************************************
       460 [3425]: ******************************************************************
       461 [3407]: ******************************************************************
       462 [3432]: *******************************************************************
       463 [3375]: *****************************************************************
       464 [3340]: *****************************************************************
       465 [3201]: **************************************************************
       466 [2562]: **************************************************
       467 [1830]: ***********************************
       468 [1158]: **********************
       469 [482 ]: *********
       470 [202 ]: ***
       471 [131 ]: **
       472 [101 ]: *
       473 [79  ]: *
       474 [66  ]: *
       475 [59  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.041552 GB

       Integrated   413 (sum) +   410 (prf) /  463 reflections on image 452
       Integrated   556 (sum) +   553 (prf) /  637 reflections on image 453
       Integrated   599 (sum) +   596 (prf) /  679 reflections on image 454
       Integrated   566 (sum) +   563 (prf) /  655 reflections on image 455
       Integrated   602 (sum) +   600 (prf) /  671 reflections on image 456
       Integrated   548 (sum) +   547 (prf) /  624 reflections on image 457
       Integrated   651 (sum) +   626 (prf) /  742 reflections on image 458
       Integrated   621 (sum) +   614 (prf) /  703 reflections on image 459
       Integrated   614 (sum) +   608 (prf) /  704 reflections on image 460
       Integrated   577 (sum) +   567 (prf) /  677 reflections on image 461
       Integrated   606 (sum) +   601 (prf) /  699 reflections on image 462
       Integrated   608 (sum) +   599 (prf) /  687 reflections on image 463
       Integrated   613 (sum) +   605 (prf) /  702 reflections on image 464
       Integrated   630 (sum) +   627 (prf) /  708 reflections on image 465
       Integrated   649 (sum) +   644 (prf) /  732 reflections on image 466
       Integrated   590 (sum) +   587 (prf) /  672 reflections on image 467
       Integrated   597 (sum) +   589 (prf) /  676 reflections on image 468
       Integrated   245 (sum) +   241 (prf) /  280 reflections on image 469
       Integrated    64 (sum) +    64 (prf) /   71 reflections on image 470
       Integrated    28 (sum) +    28 (prf) /   30 reflections on image 471
       Integrated    18 (sum) +    17 (prf) /   22 reflections on image 472
       Integrated    10 (sum) +    10 (prf) /   13 reflections on image 473
       Integrated     6 (sum) +     6 (prf) /    7 reflections on image 474
       Integrated    49 (sum) +    49 (prf) /   59 reflections on image 475
       Beginning integration job 27

       Frames: 459 -> 493

       Number of reflections
        Partial:     158
        Full:        11733
        In ice ring: 0
        Integrate:   11891
        Total:       11891

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       459 [6   ]:
       460 [13  ]:
       461 [21  ]:
       462 [38  ]:
       463 [67  ]: *
       464 [115 ]: **
       465 [248 ]: ****
       466 [847 ]: ****************
       467 [1547]: *****************************
       468 [2264]: *******************************************
       469 [2957]: *********************************************************
       470 [3207]: *************************************************************
       471 [3312]: ***************************************************************
       472 [3338]: ****************************************************************
       473 [3379]: *****************************************************************
       474 [3401]: *****************************************************************
       475 [3414]: *****************************************************************
       476 [3448]: ******************************************************************
       477 [3470]: *******************************************************************
       478 [3424]: ******************************************************************
       479 [3384]: *****************************************************************
       480 [3355]: ****************************************************************
       481 [3388]: *****************************************************************
       482 [3203]: *************************************************************
       483 [2558]: *************************************************
       484 [1892]: ************************************
       485 [1173]: **********************
       486 [472 ]: *********
       487 [226 ]: ****
       488 [168 ]: ***
       489 [126 ]: **
       490 [104 ]: **
       491 [83  ]: *
       492 [63  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0424323 GB

       Integrated   365 (sum) +   365 (prf) /  414 reflections on image 469
       Integrated   540 (sum) +   537 (prf) /  603 reflections on image 470
       Integrated   612 (sum) +   607 (prf) /  685 reflections on image 471
       Integrated   592 (sum) +   586 (prf) /  667 reflections on image 472
       Integrated   596 (sum) +   588 (prf) /  674 reflections on image 473
       Integrated   615 (sum) +   612 (prf) /  695 reflections on image 474
       Integrated   639 (sum) +   604 (prf) /  739 reflections on image 475
       Integrated   565 (sum) +   563 (prf) /  645 reflections on image 476
       Integrated   642 (sum) +   630 (prf) /  725 reflections on image 477
       Integrated   629 (sum) +   615 (prf) /  729 reflections on image 478
       Integrated   596 (sum) +   593 (prf) /  683 reflections on image 479
       Integrated   546 (sum) +   540 (prf) /  640 reflections on image 480
       Integrated   635 (sum) +   631 (prf) /  722 reflections on image 481
       Integrated   621 (sum) +   611 (prf) /  712 reflections on image 482
       Integrated   599 (sum) +   594 (prf) /  666 reflections on image 483
       Integrated   632 (sum) +   628 (prf) /  719 reflections on image 484
       Integrated   619 (sum) +   617 (prf) /  701 reflections on image 485
       Integrated   213 (sum) +   213 (prf) /  246 reflections on image 486
       Integrated    51 (sum) +    51 (prf) /   58 reflections on image 487
       Integrated    35 (sum) +    35 (prf) /   42 reflections on image 488
       Integrated    20 (sum) +    19 (prf) /   22 reflections on image 489
       Integrated    17 (sum) +    17 (prf) /   21 reflections on image 490
       Integrated    12 (sum) +    12 (prf) /   20 reflections on image 491
       Integrated    50 (sum) +    50 (prf) /   63 reflections on image 492
       Beginning integration job 28

       Frames: 476 -> 510

       Number of reflections
        Partial:     156
        Full:        11701
        In ice ring: 0
        Integrate:   11857
        Total:       11857

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       476 [5   ]:
       477 [14  ]:
       478 [23  ]:
       479 [42  ]:
       480 [65  ]: *
       481 [128 ]: **
       482 [256 ]: *****
       483 [854 ]: ****************
       484 [1574]: ******************************
       485 [2230]: *******************************************
       486 [2941]: *********************************************************
       487 [3218]: **************************************************************
       488 [3304]: ****************************************************************
       489 [3292]: ****************************************************************
       490 [3309]: ****************************************************************
       491 [3334]: *****************************************************************
       492 [3301]: ****************************************************************
       493 [3404]: ******************************************************************
       494 [3425]: *******************************************************************
       495 [3411]: ******************************************************************
       496 [3376]: ******************************************************************
       497 [3360]: *****************************************************************
       498 [3314]: ****************************************************************
       499 [3132]: *************************************************************
       500 [2549]: *************************************************
       501 [1877]: ************************************
       502 [1211]: ***********************
       503 [499 ]: *********
       504 [235 ]: ****
       505 [156 ]: ***
       506 [123 ]: **
       507 [102 ]: *
       508 [82  ]: *
       509 [70  ]: *

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0417867 GB

       Integrated   362 (sum) +   362 (prf) /  419 reflections on image 486
       Integrated   539 (sum) +   536 (prf) /  607 reflections on image 487
       Integrated   616 (sum) +   611 (prf) /  676 reflections on image 488
       Integrated   572 (sum) +   565 (prf) /  646 reflections on image 489
       Integrated   625 (sum) +   621 (prf) /  708 reflections on image 490
       Integrated   620 (sum) +   614 (prf) /  700 reflections on image 491
       Integrated   630 (sum) +   600 (prf) /  710 reflections on image 492
       Integrated   589 (sum) +   580 (prf) /  673 reflections on image 493
       Integrated   582 (sum) +   573 (prf) /  683 reflections on image 494
       Integrated   608 (sum) +   598 (prf) /  690 reflections on image 495
       Integrated   593 (sum) +   586 (prf) /  671 reflections on image 496
       Integrated   619 (sum) +   612 (prf) /  720 reflections on image 497
       Integrated   668 (sum) +   664 (prf) /  744 reflections on image 498
       Integrated   585 (sum) +   580 (prf) /  661 reflections on image 499
       Integrated   580 (sum) +   574 (prf) /  672 reflections on image 500
       Integrated   591 (sum) +   584 (prf) /  666 reflections on image 501
       Integrated   622 (sum) +   617 (prf) /  712 reflections on image 502
       Integrated   230 (sum) +   228 (prf) /  264 reflections on image 503
       Integrated    69 (sum) +    68 (prf) /   79 reflections on image 504
       Integrated    29 (sum) +    29 (prf) /   33 reflections on image 505
       Integrated    17 (sum) +    17 (prf) /   21 reflections on image 506
       Integrated    20 (sum) +    20 (prf) /   20 reflections on image 507
       Integrated    12 (sum) +    12 (prf) /   12 reflections on image 508
       Integrated    61 (sum) +    61 (prf) /   70 reflections on image 509
       Beginning integration job 29

       Frames: 493 -> 527

       Number of reflections
        Partial:     88
        Full:        11024
        In ice ring: 0
        Integrate:   11112
        Total:       11112

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       493 [8   ]:
       494 [19  ]:
       495 [35  ]:
       496 [48  ]:
       497 [71  ]: *
       498 [125 ]: **
       499 [246 ]: ****
       500 [847 ]: ****************
       501 [1545]: ******************************
       502 [2257]: ********************************************
       503 [2942]: *********************************************************
       504 [3256]: ****************************************************************
       505 [3335]: *****************************************************************
       506 [3319]: *****************************************************************
       507 [3329]: *****************************************************************
       508 [3316]: *****************************************************************
       509 [3405]: *******************************************************************
       510 [3377]: ******************************************************************
       511 [3405]: *******************************************************************
       512 [3396]: ******************************************************************
       513 [3326]: *****************************************************************
       514 [3242]: ***************************************************************
       515 [3142]: *************************************************************
       516 [2516]: *************************************************
       517 [1801]: ***********************************
       518 [1145]: **********************
       519 [445 ]: ********
       520 [173 ]: ***
       521 [89  ]: *
       522 [55  ]: *
       523 [38  ]:
       524 [21  ]:
       525 [16  ]:
       526 [6   ]:

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0408313 GB

       Integrated   347 (sum) +   347 (prf) /  402 reflections on image 503
       Integrated   525 (sum) +   521 (prf) /  607 reflections on image 504
       Integrated   618 (sum) +   616 (prf) /  696 reflections on image 505
       Integrated   577 (sum) +   576 (prf) /  661 reflections on image 506
       Integrated   610 (sum) +   605 (prf) /  685 reflections on image 507
       Integrated   569 (sum) +   568 (prf) /  647 reflections on image 508
       Integrated   669 (sum) +   658 (prf) /  752 reflections on image 509
       Integrated   572 (sum) +   560 (prf) /  652 reflections on image 510
       Integrated   583 (sum) +   570 (prf) /  673 reflections on image 511
       Integrated   639 (sum) +   632 (prf) /  740 reflections on image 512
       Integrated   631 (sum) +   625 (prf) /  721 reflections on image 513
       Integrated   580 (sum) +   573 (prf) /  654 reflections on image 514
       Integrated   626 (sum) +   618 (prf) /  706 reflections on image 515
       Integrated   638 (sum) +   633 (prf) /  715 reflections on image 516
       Integrated   582 (sum) +   577 (prf) /  656 reflections on image 517
       Integrated   622 (sum) +   615 (prf) /  700 reflections on image 518
       Integrated   237 (sum) +   234 (prf) /  272 reflections on image 519
       Integrated    71 (sum) +    71 (prf) /   84 reflections on image 520
       Integrated    30 (sum) +    30 (prf) /   34 reflections on image 521
       Integrated    14 (sum) +    14 (prf) /   17 reflections on image 522
       Integrated    15 (sum) +    15 (prf) /   17 reflections on image 523
       Integrated     4 (sum) +     4 (prf) /    5 reflections on image 524
       Integrated     8 (sum) +     8 (prf) /   10 reflections on image 525
       Integrated     4 (sum) +     4 (prf) /    6 reflections on image 526
       Beginning integration job 30

       Frames: 510 -> 540

       Number of reflections
        Partial:     1333
        Full:        14230
        In ice ring: 0
        Integrate:   15563
        Total:       15563

       The following histogram shows the number of reflections predicted
       to have all or part of their intensity on each frame.

       510 [52  ]: *
       511 [67  ]: *
       512 [80  ]: *
       513 [103 ]: *
       514 [151 ]: **
       515 [307 ]: *****
       516 [925 ]: *****************
       517 [1642]: *******************************
       518 [2292]: ********************************************
       519 [3022]: **********************************************************
       520 [3289]: ***************************************************************
       521 [3306]: ***************************************************************
       522 [3329]: ****************************************************************
       523 [3362]: *****************************************************************
       524 [3428]: ******************************************************************
       525 [3464]: *******************************************************************
       526 [3430]: ******************************************************************
       527 [3401]: *****************************************************************
       528 [3444]: ******************************************************************
       529 [3372]: *****************************************************************
       530 [3381]: *****************************************************************
       531 [3375]: *****************************************************************
       532 [3316]: ****************************************************************
       533 [3333]: ****************************************************************
       534 [3417]: ******************************************************************
       535 [3411]: *****************************************************************
       536 [3395]: *****************************************************************
       537 [3207]: **************************************************************
       538 [2737]: ****************************************************
       539 [2025]: ***************************************

       Memory usage:
        Total system memory: 16.726 GB
        Limit shoebox memory: 3.13612 GB
        Required shoebox memory: 0.0422526 GB

       Integrated   378 (sum) +   377 (prf) /  431 reflections on image 519
       Integrated   546 (sum) +   543 (prf) /  639 reflections on image 520
       Integrated   590 (sum) +   587 (prf) /  663 reflections on image 521
       Integrated   609 (sum) +   606 (prf) /  687 reflections on image 522
       Integrated   572 (sum) +   565 (prf) /  655 reflections on image 523
       Integrated   580 (sum) +   577 (prf) /  658 reflections on image 524
       Integrated   599 (sum) +   595 (prf) /  697 reflections on image 525
       Integrated   656 (sum) +   649 (prf) /  755 reflections on image 526
       Integrated   579 (sum) +   570 (prf) /  673 reflections on image 527
       Integrated   641 (sum) +   636 (prf) /  731 reflections on image 528
       Integrated   600 (sum) +   594 (prf) /  679 reflections on image 529
       Integrated   610 (sum) +   601 (prf) /  702 reflections on image 530
       Integrated   647 (sum) +   639 (prf) /  724 reflections on image 531
       Integrated   615 (sum) +   612 (prf) /  701 reflections on image 532
       Integrated   571 (sum) +   567 (prf) /  631 reflections on image 533
       Integrated   623 (sum) +   619 (prf) /  693 reflections on image 534
       Integrated   590 (sum) +   585 (prf) /  663 reflections on image 535
       Integrated   654 (sum) +   648 (prf) /  731 reflections on image 536
       Integrated   620 (sum) +   614 (prf) /  713 reflections on image 537
       Integrated   630 (sum) +   624 (prf) /  712 reflections on image 538
       Integrated  1774 (sum) +  1758 (prf) / 2025 reflections on image 539

       Summary vs image number
       ------------------------------------------------------------------------------------------------------
       ID | Image | # full | # part | # over | # ice | # sum | # prf | <Ibg> | <I/sigI> | <I/sigI> | <CC prf>
          |       |        |        |        |       |       |       |       |  (sum)   |  (prf)   |
       ------------------------------------------------------------------------------------------------------
       0  | 0     | 0      | 706    | 0      | 0     | 637   | 633   | 0.15  | 5.17     | 5.39     | 0.45
       0  | 1     | 292    | 413    | 0      | 0     | 623   | 621   | 0.15  | 5.80     | 6.04     | 0.46
       0  | 2     | 664    | 81     | 0      | 0     | 660   | 657   | 0.15  | 5.26     | 5.62     | 0.46
       0  | 3     | 646    | 42     | 0      | 0     | 611   | 607   | 0.15  | 5.25     | 5.56     | 0.45
       0  | 4     | 657    | 23     | 0      | 0     | 604   | 596   | 0.15  | 6.20     | 6.55     | 0.46
       0  | 5     | 629    | 20     | 0      | 0     | 555   | 549   | 0.15  | 5.58     | 5.92     | 0.46
       0  | 6     | 700    | 7      | 0      | 0     | 616   | 610   | 0.15  | 5.82     | 6.15     | 0.44
       0  | 7     | 667    | 6      | 0      | 0     | 588   | 583   | 0.16  | 6.12     | 6.41     | 0.46
       0  | 8     | 695    | 11     | 0      | 0     | 603   | 597   | 0.15  | 5.31     | 5.66     | 0.44
       0  | 9     | 680    | 10     | 0      | 0     | 590   | 585   | 0.15  | 5.29     | 5.63     | 0.46
       0  | 10    | 696    | 6      | 0      | 0     | 619   | 613   | 0.16  | 5.55     | 5.91     | 0.47
       0  | 11    | 699    | 1      | 0      | 0     | 600   | 593   | 0.14  | 4.82     | 5.19     | 0.45
       0  | 12    | 634    | 4      | 0      | 0     | 545   | 544   | 0.15  | 5.59     | 5.87     | 0.45
       0  | 13    | 683    | 3      | 0      | 0     | 614   | 608   | 0.16  | 5.64     | 5.98     | 0.45
       0  | 14    | 653    | 3      | 0      | 0     | 563   | 562   | 0.16  | 5.68     | 5.95     | 0.46
       0  | 15    | 649    | 10     | 0      | 0     | 580   | 575   | 0.16  | 6.52     | 6.84     | 0.46
       0  | 16    | 725    | 5      | 0      | 0     | 631   | 625   | 0.15  | 5.27     | 5.59     | 0.45
       0  | 17    | 701    | 2      | 0      | 0     | 624   | 620   | 0.16  | 5.73     | 6.08     | 0.48
       0  | 18    | 685    | 8      | 0      | 0     | 620   | 611   | 0.15  | 5.11     | 5.51     | 0.48
       0  | 19    | 760    | 7      | 0      | 0     | 681   | 672   | 0.15  | 5.53     | 5.89     | 0.47
       0  | 20    | 708    | 4      | 0      | 0     | 622   | 615   | 0.16  | 6.39     | 6.72     | 0.46
       0  | 21    | 697    | 4      | 0      | 0     | 641   | 634   | 0.15  | 5.54     | 5.86     | 0.45
       0  | 22    | 656    | 13     | 0      | 0     | 595   | 588   | 0.16  | 6.16     | 6.53     | 0.46
       0  | 23    | 707    | 6      | 0      | 0     | 631   | 625   | 0.16  | 6.06     | 6.43     | 0.46
       0  | 24    | 691    | 11     | 0      | 0     | 619   | 610   | 0.15  | 5.31     | 5.68     | 0.45
       0  | 25    | 664    | 19     | 0      | 0     | 586   | 581   | 0.15  | 5.48     | 5.80     | 0.46
       0  | 26    | 706    | 14     | 0      | 0     | 629   | 621   | 0.15  | 5.74     | 6.04     | 0.46
       0  | 27    | 725    | 13     | 0      | 0     | 639   | 636   | 0.16  | 5.95     | 6.23     | 0.46
       0  | 28    | 693    | 9      | 0      | 0     | 609   | 602   | 0.16  | 5.85     | 6.18     | 0.46
       0  | 29    | 718    | 10     | 0      | 0     | 633   | 630   | 0.16  | 6.04     | 6.35     | 0.47
       0  | 30    | 725    | 4      | 0      | 0     | 656   | 655   | 0.15  | 5.63     | 5.93     | 0.44
       0  | 31    | 654    | 4      | 0      | 0     | 567   | 566   | 0.16  | 5.98     | 6.28     | 0.47
       0  | 32    | 655    | 5      | 0      | 0     | 572   | 565   | 0.15  | 5.04     | 5.35     | 0.44
       0  | 33    | 680    | 6      | 0      | 0     | 592   | 588   | 0.15  | 5.65     | 5.94     | 0.43
       0  | 34    | 679    | 10     | 0      | 0     | 594   | 587   | 0.15  | 5.50     | 5.85     | 0.46
       0  | 35    | 713    | 3      | 0      | 0     | 633   | 623   | 0.15  | 4.94     | 5.32     | 0.43
       0  | 36    | 661    | 3      | 0      | 0     | 598   | 593   | 0.15  | 5.49     | 5.83     | 0.45
       0  | 37    | 685    | 6      | 0      | 0     | 613   | 609   | 0.16  | 5.75     | 6.04     | 0.46
       0  | 38    | 688    | 15     | 0      | 0     | 616   | 610   | 0.16  | 5.86     | 6.16     | 0.44
       0  | 39    | 712    | 19     | 0      | 0     | 645   | 635   | 0.15  | 5.17     | 5.52     | 0.44
       0  | 40    | 662    | 9      | 0      | 0     | 593   | 589   | 0.15  | 5.57     | 5.85     | 0.43
       0  | 41    | 622    | 12     | 0      | 0     | 547   | 541   | 0.16  | 5.57     | 5.90     | 0.44
       0  | 42    | 666    | 21     | 0      | 0     | 603   | 593   | 0.16  | 5.65     | 6.05     | 0.45
       0  | 43    | 679    | 25     | 0      | 0     | 617   | 608   | 0.16  | 6.45     | 6.77     | 0.47
       0  | 44    | 725    | 18     | 0      | 0     | 660   | 657   | 0.16  | 5.79     | 6.13     | 0.46
       0  | 45    | 730    | 13     | 0      | 0     | 648   | 643   | 0.16  | 5.94     | 6.28     | 0.47
       0  | 46    | 703    | 6      | 0      | 0     | 612   | 607   | 0.16  | 5.71     | 5.99     | 0.46
       0  | 47    | 669    | 11     | 0      | 0     | 562   | 558   | 0.15  | 5.21     | 5.52     | 0.46
       0  | 48    | 720    | 10     | 0      | 0     | 647   | 638   | 0.16  | 5.79     | 6.11     | 0.45
       0  | 49    | 729    | 7      | 0      | 0     | 644   | 634   | 0.15  | 5.70     | 6.03     | 0.46
       0  | 50    | 738    | 3      | 0      | 0     | 648   | 644   | 0.16  | 5.60     | 5.93     | 0.46
       0  | 51    | 651    | 5      | 0      | 0     | 580   | 579   | 0.15  | 5.91     | 6.22     | 0.46
       0  | 52    | 664    | 4      | 0      | 0     | 585   | 579   | 0.15  | 5.57     | 5.87     | 0.44
       0  | 53    | 701    | 13     | 0      | 0     | 625   | 616   | 0.15  | 5.74     | 6.11     | 0.44
       0  | 54    | 681    | 13     | 0      | 0     | 611   | 599   | 0.14  | 4.99     | 5.32     | 0.43
       0  | 55    | 639    | 6      | 0      | 0     | 569   | 566   | 0.16  | 5.64     | 5.92     | 0.46
       0  | 56    | 690    | 20     | 0      | 0     | 612   | 604   | 0.15  | 4.95     | 5.32     | 0.44
       0  | 57    | 721    | 27     | 0      | 0     | 668   | 657   | 0.15  | 6.30     | 6.57     | 0.44
       0  | 58    | 684    | 15     | 0      | 0     | 607   | 602   | 0.15  | 5.28     | 5.59     | 0.44
       0  | 59    | 683    | 35     | 0      | 0     | 622   | 610   | 0.17  | 6.31     | 6.65     | 0.48
       0  | 60    | 713    | 9      | 0      | 0     | 640   | 633   | 0.16  | 5.71     | 6.04     | 0.44
       0  | 61    | 658    | 7      | 0      | 0     | 575   | 572   | 0.15  | 5.30     | 5.56     | 0.43
       0  | 62    | 681    | 3      | 0      | 0     | 595   | 592   | 0.16  | 6.17     | 6.46     | 0.45
       0  | 63    | 640    | 22     | 0      | 0     | 569   | 567   | 0.15  | 5.54     | 5.84     | 0.45
       0  | 64    | 670    | 3      | 0      | 0     | 598   | 592   | 0.15  | 6.01     | 6.28     | 0.45
       0  | 65    | 692    | 0      | 0      | 0     | 611   | 606   | 0.16  | 6.36     | 6.66     | 0.48
       0  | 66    | 710    | 16     | 0      | 0     | 637   | 634   | 0.16  | 6.47     | 6.73     | 0.44
       0  | 67    | 705    | 4      | 0      | 0     | 632   | 626   | 0.16  | 5.80     | 6.12     | 0.45
       0  | 68    | 708    | 4      | 0      | 0     | 628   | 620   | 0.16  | 5.61     | 5.95     | 0.45
       0  | 69    | 712    | 9      | 0      | 0     | 634   | 629   | 0.15  | 5.41     | 5.77     | 0.45
       0  | 70    | 711    | 19     | 0      | 0     | 619   | 612   | 0.16  | 6.29     | 6.56     | 0.47
       0  | 71    | 705    | 7      | 0      | 0     | 624   | 619   | 0.15  | 5.32     | 5.62     | 0.45
       0  | 72    | 723    | 6      | 0      | 0     | 624   | 617   | 0.15  | 5.18     | 5.53     | 0.45
       0  | 73    | 692    | 3      | 0      | 0     | 605   | 599   | 0.17  | 5.97     | 6.35     | 0.48
       0  | 74    | 680    | 26     | 0      | 0     | 606   | 600   | 0.16  | 6.10     | 6.37     | 0.45
       0  | 75    | 672    | 9      | 0      | 0     | 594   | 587   | 0.15  | 5.24     | 5.61     | 0.46
       0  | 76    | 686    | 43     | 0      | 0     | 618   | 594   | 0.15  | 5.78     | 6.21     | 0.44
       0  | 77    | 681    | 15     | 0      | 0     | 605   | 598   | 0.15  | 5.19     | 5.54     | 0.44
       0  | 78    | 645    | 12     | 0      | 0     | 569   | 564   | 0.16  | 6.14     | 6.43     | 0.45
       0  | 79    | 656    | 3      | 0      | 0     | 587   | 582   | 0.15  | 5.57     | 5.89     | 0.45
       0  | 80    | 685    | 16     | 0      | 0     | 607   | 601   | 0.15  | 5.92     | 6.22     | 0.45
       0  | 81    | 712    | 6      | 0      | 0     | 618   | 612   | 0.15  | 5.26     | 5.57     | 0.43
       0  | 82    | 653    | 14     | 0      | 0     | 585   | 577   | 0.16  | 6.19     | 6.51     | 0.46
       0  | 83    | 672    | 3      | 0      | 0     | 591   | 583   | 0.16  | 5.77     | 6.15     | 0.47
       0  | 84    | 704    | 4      | 0      | 0     | 640   | 634   | 0.16  | 5.71     | 6.01     | 0.46
       0  | 85    | 758    | 4      | 0      | 0     | 668   | 664   | 0.16  | 5.83     | 6.17     | 0.45
       0  | 86    | 740    | 12     | 0      | 0     | 659   | 649   | 0.16  | 5.62     | 5.96     | 0.46
       0  | 87    | 682    | 7      | 0      | 0     | 621   | 615   | 0.15  | 5.69     | 5.99     | 0.46
       0  | 88    | 686    | 9      | 0      | 0     | 608   | 601   | 0.15  | 5.57     | 5.90     | 0.45
       0  | 89    | 702    | 3      | 0      | 0     | 616   | 613   | 0.15  | 5.23     | 5.53     | 0.46
       0  | 90    | 667    | 6      | 0      | 0     | 596   | 591   | 0.14  | 4.99     | 5.28     | 0.43
       0  | 91    | 717    | 12     | 0      | 0     | 645   | 635   | 0.16  | 5.53     | 5.84     | 0.43
       0  | 92    | 671    | 30     | 0      | 0     | 620   | 608   | 0.15  | 5.35     | 5.68     | 0.42
       0  | 93    | 616    | 9      | 0      | 0     | 542   | 534   | 0.15  | 5.30     | 5.66     | 0.46
       0  | 94    | 703    | 14     | 0      | 0     | 621   | 617   | 0.16  | 6.15     | 6.44     | 0.46
       0  | 95    | 712    | 7      | 0      | 0     | 639   | 635   | 0.16  | 5.85     | 6.11     | 0.45
       0  | 96    | 691    | 12     | 0      | 0     | 608   | 604   | 0.16  | 6.17     | 6.40     | 0.44
       0  | 97    | 651    | 29     | 0      | 0     | 592   | 585   | 0.16  | 6.11     | 6.37     | 0.45
       0  | 98    | 715    | 18     | 0      | 0     | 638   | 634   | 0.16  | 6.45     | 6.76     | 0.47
       0  | 99    | 648    | 7      | 0      | 0     | 574   | 563   | 0.17  | 5.89     | 6.28     | 0.46
       0  | 100   | 691    | 19     | 0      | 0     | 605   | 600   | 0.16  | 5.91     | 6.21     | 0.45
       0  | 101   | 755    | 6      | 0      | 0     | 657   | 652   | 0.16  | 5.67     | 6.00     | 0.46
       0  | 102   | 674    | 7      | 0      | 0     | 602   | 595   | 0.15  | 4.66     | 5.01     | 0.42
       0  | 103   | 694    | 9      | 0      | 0     | 604   | 596   | 0.16  | 5.54     | 5.91     | 0.46
       0  | 104   | 696    | 20     | 0      | 0     | 612   | 606   | 0.16  | 6.06     | 6.39     | 0.46
       0  | 105   | 701    | 11     | 0      | 0     | 602   | 598   | 0.16  | 5.61     | 5.92     | 0.46
       0  | 106   | 693    | 7      | 0      | 0     | 609   | 604   | 0.16  | 5.33     | 5.63     | 0.44
       0  | 107   | 691    | 25     | 0      | 0     | 623   | 614   | 0.16  | 6.51     | 6.81     | 0.45
       0  | 108   | 688    | 6      | 0      | 0     | 614   | 607   | 0.15  | 4.76     | 5.09     | 0.42
       0  | 109   | 676    | 15     | 0      | 0     | 601   | 590   | 0.15  | 5.48     | 5.83     | 0.45
       0  | 110   | 701    | 29     | 0      | 0     | 636   | 623   | 0.15  | 5.58     | 5.95     | 0.44
       0  | 111   | 685    | 11     | 0      | 0     | 612   | 606   | 0.16  | 5.69     | 6.03     | 0.46
       0  | 112   | 683    | 15     | 0      | 0     | 626   | 621   | 0.16  | 5.83     | 6.12     | 0.44
       0  | 113   | 625    | 26     | 0      | 0     | 559   | 555   | 0.16  | 6.20     | 6.53     | 0.46
       0  | 114   | 746    | 7      | 0      | 0     | 664   | 658   | 0.17  | 6.57     | 6.88     | 0.46
       0  | 115   | 647    | 9      | 0      | 0     | 570   | 566   | 0.16  | 5.34     | 5.65     | 0.45
       0  | 116   | 682    | 12     | 0      | 0     | 609   | 599   | 0.16  | 5.83     | 6.20     | 0.47
       0  | 117   | 744    | 11     | 0      | 0     | 661   | 654   | 0.17  | 6.41     | 6.71     | 0.47
       0  | 118   | 722    | 4      | 0      | 0     | 623   | 615   | 0.15  | 5.32     | 5.63     | 0.46
       0  | 119   | 681    | 8      | 0      | 0     | 598   | 595   | 0.16  | 5.36     | 5.69     | 0.47
       0  | 120   | 686    | 7      | 0      | 0     | 590   | 589   | 0.15  | 5.20     | 5.45     | 0.43
       0  | 121   | 702    | 8      | 0      | 0     | 631   | 627   | 0.15  | 5.97     | 6.31     | 0.46
       0  | 122   | 711    | 19     | 0      | 0     | 636   | 628   | 0.15  | 5.59     | 5.86     | 0.44
       0  | 123   | 658    | 3      | 0      | 0     | 588   | 578   | 0.16  | 5.26     | 5.62     | 0.44
       0  | 124   | 679    | 17     | 0      | 0     | 624   | 617   | 0.15  | 6.01     | 6.32     | 0.45
       0  | 125   | 692    | 9      | 0      | 0     | 605   | 600   | 0.16  | 5.30     | 5.56     | 0.44
       0  | 126   | 650    | 12     | 0      | 0     | 584   | 575   | 0.16  | 5.44     | 5.79     | 0.46
       0  | 127   | 698    | 21     | 0      | 0     | 634   | 626   | 0.16  | 5.86     | 6.24     | 0.48
       0  | 128   | 683    | 0      | 0      | 0     | 611   | 606   | 0.16  | 5.30     | 5.67     | 0.47
       0  | 129   | 603    | 15     | 0      | 0     | 535   | 526   | 0.16  | 5.99     | 6.31     | 0.45
       0  | 130   | 725    | 21     | 0      | 0     | 656   | 643   | 0.16  | 5.82     | 6.17     | 0.45
       0  | 131   | 709    | 3      | 0      | 0     | 623   | 621   | 0.15  | 5.40     | 5.70     | 0.44
       0  | 132   | 702    | 3      | 0      | 0     | 602   | 597   | 0.16  | 5.53     | 5.83     | 0.45
       0  | 133   | 736    | 15     | 0      | 0     | 654   | 646   | 0.16  | 5.70     | 6.06     | 0.47
       0  | 134   | 731    | 3      | 0      | 0     | 634   | 629   | 0.16  | 5.86     | 6.16     | 0.46
       0  | 135   | 669    | 0      | 0      | 0     | 584   | 581   | 0.15  | 5.77     | 6.06     | 0.44
       0  | 136   | 710    | 4      | 0      | 0     | 624   | 614   | 0.16  | 5.58     | 5.95     | 0.44
       0  | 137   | 692    | 7      | 0      | 0     | 619   | 614   | 0.16  | 5.78     | 6.13     | 0.45
       0  | 138   | 699    | 31     | 0      | 0     | 641   | 633   | 0.15  | 5.39     | 5.70     | 0.45
       0  | 139   | 680    | 3      | 0      | 0     | 591   | 586   | 0.17  | 6.26     | 6.63     | 0.47
       0  | 140   | 653    | 12     | 0      | 0     | 596   | 588   | 0.16  | 5.36     | 5.70     | 0.45
       0  | 141   | 663    | 22     | 0      | 0     | 582   | 579   | 0.16  | 5.46     | 5.76     | 0.46
       0  | 142   | 658    | 10     | 0      | 0     | 584   | 576   | 0.15  | 5.59     | 5.92     | 0.44
       0  | 143   | 620    | 6      | 0      | 0     | 561   | 555   | 0.16  | 5.48     | 5.84     | 0.47
       0  | 144   | 729    | 21     | 0      | 0     | 653   | 643   | 0.16  | 5.39     | 5.67     | 0.42
       0  | 145   | 738    | 35     | 0      | 0     | 661   | 648   | 0.15  | 4.64     | 5.05     | 0.43
       0  | 146   | 707    | 16     | 0      | 0     | 642   | 633   | 0.15  | 5.38     | 5.70     | 0.42
       0  | 147   | 735    | 10     | 0      | 0     | 643   | 637   | 0.17  | 5.82     | 6.14     | 0.46
       0  | 148   | 665    | 9      | 0      | 0     | 598   | 590   | 0.15  | 5.09     | 5.40     | 0.43
       0  | 149   | 674    | 10     | 0      | 0     | 591   | 586   | 0.16  | 6.27     | 6.59     | 0.46
       0  | 150   | 708    | 13     | 0      | 0     | 645   | 638   | 0.16  | 6.11     | 6.45     | 0.45
       0  | 151   | 676    | 21     | 0      | 0     | 602   | 596   | 0.17  | 6.22     | 6.57     | 0.47
       0  | 152   | 697    | 30     | 0      | 0     | 604   | 596   | 0.16  | 5.68     | 6.05     | 0.48
       0  | 153   | 690    | 4      | 0      | 0     | 615   | 612   | 0.16  | 6.18     | 6.45     | 0.45
       0  | 154   | 640    | 12     | 0      | 0     | 588   | 582   | 0.17  | 6.16     | 6.50     | 0.47
       0  | 155   | 669    | 3      | 0      | 0     | 605   | 601   | 0.16  | 5.82     | 6.12     | 0.46
       0  | 156   | 698    | 16     | 0      | 0     | 614   | 606   | 0.16  | 5.60     | 5.92     | 0.44
       0  | 157   | 721    | 3      | 0      | 0     | 628   | 622   | 0.16  | 5.76     | 6.09     | 0.45
       0  | 158   | 729    | 3      | 0      | 0     | 644   | 637   | 0.16  | 5.97     | 6.29     | 0.47
       0  | 159   | 704    | 9      | 0      | 0     | 621   | 612   | 0.16  | 5.03     | 5.37     | 0.44
       0  | 160   | 713    | 24     | 0      | 0     | 642   | 628   | 0.16  | 6.02     | 6.39     | 0.47
       0  | 161   | 664    | 25     | 0      | 0     | 609   | 601   | 0.16  | 5.90     | 6.26     | 0.45
       0  | 162   | 697    | 32     | 0      | 0     | 638   | 620   | 0.16  | 5.80     | 6.21     | 0.47
       0  | 163   | 687    | 32     | 0      | 0     | 617   | 611   | 0.16  | 5.79     | 6.09     | 0.45
       0  | 164   | 688    | 49     | 0      | 0     | 622   | 613   | 0.17  | 6.35     | 6.64     | 0.46
       0  | 165   | 666    | 18     | 0      | 0     | 593   | 584   | 0.17  | 5.74     | 6.08     | 0.46
       0  | 166   | 686    | 15     | 0      | 0     | 601   | 596   | 0.16  | 5.18     | 5.46     | 0.44
       0  | 167   | 642    | 3      | 0      | 0     | 562   | 559   | 0.17  | 6.22     | 6.54     | 0.49
       0  | 168   | 683    | 4      | 0      | 0     | 601   | 592   | 0.16  | 6.03     | 6.36     | 0.46
       0  | 169   | 732    | 4      | 0      | 0     | 643   | 638   | 0.15  | 5.42     | 5.70     | 0.44
       0  | 170   | 660    | 4      | 0      | 0     | 585   | 581   | 0.16  | 5.56     | 5.90     | 0.44
       0  | 171   | 732    | 10     | 0      | 0     | 653   | 647   | 0.16  | 5.68     | 5.98     | 0.45
       0  | 172   | 715    | 7      | 0      | 0     | 643   | 637   | 0.16  | 6.16     | 6.39     | 0.46
       0  | 173   | 685    | 7      | 0      | 0     | 618   | 613   | 0.17  | 6.03     | 6.34     | 0.46
       0  | 174   | 686    | 9      | 0      | 0     | 617   | 613   | 0.16  | 5.68     | 5.97     | 0.44
       0  | 175   | 667    | 13     | 0      | 0     | 600   | 593   | 0.15  | 5.18     | 5.52     | 0.45
       0  | 176   | 702    | 15     | 0      | 0     | 622   | 611   | 0.16  | 5.34     | 5.72     | 0.46
       0  | 177   | 630    | 14     | 0      | 0     | 570   | 561   | 0.16  | 6.08     | 6.43     | 0.46
       0  | 178   | 673    | 22     | 0      | 0     | 590   | 577   | 0.16  | 5.54     | 5.95     | 0.46
       0  | 179   | 667    | 19     | 0      | 0     | 594   | 588   | 0.17  | 6.21     | 6.55     | 0.47
       0  | 180   | 682    | 6      | 0      | 0     | 603   | 597   | 0.17  | 5.13     | 5.45     | 0.45
       0  | 181   | 724    | 12     | 0      | 0     | 645   | 641   | 0.17  | 6.62     | 6.82     | 0.43
       0  | 182   | 711    | 12     | 0      | 0     | 631   | 623   | 0.16  | 5.81     | 6.14     | 0.46
       0  | 183   | 704    | 13     | 0      | 0     | 624   | 617   | 0.16  | 5.33     | 5.63     | 0.44
       0  | 184   | 678    | 12     | 0      | 0     | 610   | 603   | 0.15  | 5.10     | 5.46     | 0.44
       0  | 185   | 716    | 7      | 0      | 0     | 645   | 638   | 0.16  | 5.74     | 6.09     | 0.46
       0  | 186   | 697    | 6      | 0      | 0     | 584   | 578   | 0.17  | 5.91     | 6.18     | 0.45
       0  | 187   | 676    | 0      | 0      | 0     | 593   | 588   | 0.17  | 6.03     | 6.34     | 0.46
       0  | 188   | 709    | 16     | 0      | 0     | 640   | 632   | 0.15  | 4.66     | 5.04     | 0.45
       0  | 189   | 679    | 21     | 0      | 0     | 618   | 609   | 0.17  | 5.71     | 6.08     | 0.47
       0  | 190   | 648    | 0      | 0      | 0     | 576   | 569   | 0.16  | 5.91     | 6.21     | 0.45
       0  | 191   | 678    | 7      | 0      | 0     | 603   | 597   | 0.16  | 6.08     | 6.41     | 0.46
       0  | 192   | 701    | 14     | 0      | 0     | 639   | 635   | 0.17  | 6.10     | 6.38     | 0.46
       0  | 193   | 687    | 13     | 0      | 0     | 624   | 617   | 0.16  | 5.76     | 6.09     | 0.45
       0  | 194   | 681    | 6      | 0      | 0     | 590   | 584   | 0.16  | 5.50     | 5.81     | 0.42
       0  | 195   | 710    | 25     | 0      | 0     | 629   | 618   | 0.17  | 6.15     | 6.50     | 0.47
       0  | 196   | 693    | 22     | 0      | 0     | 623   | 615   | 0.16  | 5.12     | 5.50     | 0.43
       0  | 197   | 767    | 5      | 0      | 0     | 667   | 664   | 0.16  | 5.12     | 5.46     | 0.46
       0  | 198   | 652    | 19     | 0      | 0     | 585   | 583   | 0.17  | 5.89     | 6.14     | 0.45
       0  | 199   | 683    | 13     | 0      | 0     | 615   | 608   | 0.16  | 5.65     | 5.95     | 0.44
       0  | 200   | 664    | 6      | 0      | 0     | 584   | 581   | 0.17  | 6.08     | 6.47     | 0.48
       0  | 201   | 689    | 20     | 0      | 0     | 606   | 605   | 0.16  | 5.42     | 5.71     | 0.45
       0  | 202   | 659    | 31     | 0      | 0     | 584   | 581   | 0.17  | 6.26     | 6.56     | 0.46
       0  | 203   | 701    | 0      | 0      | 0     | 604   | 597   | 0.16  | 5.94     | 6.25     | 0.44
       0  | 204   | 701    | 7      | 0      | 0     | 622   | 620   | 0.16  | 5.60     | 5.86     | 0.44
       0  | 205   | 668    | 18     | 0      | 0     | 603   | 594   | 0.17  | 5.91     | 6.27     | 0.46
       0  | 206   | 723    | 11     | 0      | 0     | 651   | 644   | 0.17  | 5.28     | 5.58     | 0.44
       0  | 207   | 707    | 7      | 0      | 0     | 621   | 618   | 0.16  | 5.35     | 5.66     | 0.46
       0  | 208   | 708    | 19     | 0      | 0     | 648   | 639   | 0.16  | 5.33     | 5.67     | 0.45
       0  | 209   | 691    | 24     | 0      | 0     | 598   | 589   | 0.17  | 6.05     | 6.39     | 0.45
       0  | 210   | 673    | 6      | 0      | 0     | 574   | 567   | 0.16  | 5.55     | 5.92     | 0.44
       0  | 211   | 693    | 14     | 0      | 0     | 618   | 613   | 0.17  | 5.88     | 6.23     | 0.47
       0  | 212   | 670    | 9      | 0      | 0     | 592   | 586   | 0.17  | 6.30     | 6.66     | 0.48
       0  | 213   | 663    | 19     | 0      | 0     | 602   | 593   | 0.17  | 5.66     | 6.00     | 0.46
       0  | 214   | 692    | 6      | 0      | 0     | 619   | 612   | 0.17  | 5.99     | 6.32     | 0.46
       0  | 215   | 680    | 14     | 0      | 0     | 600   | 591   | 0.17  | 5.76     | 6.09     | 0.44
       0  | 216   | 711    | 9      | 0      | 0     | 643   | 635   | 0.17  | 5.97     | 6.28     | 0.44
       0  | 217   | 709    | 9      | 0      | 0     | 636   | 630   | 0.17  | 5.32     | 5.65     | 0.45
       0  | 218   | 680    | 3      | 0      | 0     | 597   | 591   | 0.15  | 5.24     | 5.55     | 0.43
       0  | 219   | 718    | 3      | 0      | 0     | 633   | 628   | 0.16  | 5.30     | 5.64     | 0.44
       0  | 220   | 707    | 13     | 0      | 0     | 646   | 639   | 0.16  | 5.37     | 5.71     | 0.44
       0  | 221   | 680    | 3      | 0      | 0     | 602   | 598   | 0.17  | 5.87     | 6.25     | 0.47
       0  | 222   | 670    | 21     | 0      | 0     | 580   | 573   | 0.17  | 6.08     | 6.42     | 0.46
       0  | 223   | 686    | 14     | 0      | 0     | 610   | 605   | 0.17  | 5.63     | 5.97     | 0.47
       0  | 224   | 650    | 3      | 0      | 0     | 589   | 586   | 0.17  | 6.17     | 6.46     | 0.46
       0  | 225   | 672    | 27     | 0      | 0     | 602   | 599   | 0.17  | 5.72     | 6.04     | 0.46
       0  | 226   | 708    | 22     | 0      | 0     | 646   | 633   | 0.17  | 6.23     | 6.55     | 0.43
       0  | 227   | 731    | 9      | 0      | 0     | 640   | 635   | 0.17  | 5.91     | 6.23     | 0.47
       0  | 228   | 656    | 39     | 0      | 0     | 582   | 569   | 0.16  | 4.99     | 5.35     | 0.44
       0  | 229   | 712    | 6      | 0      | 0     | 638   | 629   | 0.16  | 4.94     | 5.35     | 0.45
       0  | 230   | 707    | 12     | 0      | 0     | 639   | 633   | 0.16  | 5.46     | 5.79     | 0.46
       0  | 231   | 714    | 17     | 0      | 0     | 622   | 616   | 0.17  | 5.44     | 5.75     | 0.44
       0  | 232   | 662    | 0      | 0      | 0     | 579   | 578   | 0.17  | 6.23     | 6.56     | 0.47
       0  | 233   | 686    | 10     | 0      | 0     | 606   | 601   | 0.16  | 5.48     | 5.90     | 0.46
       0  | 234   | 662    | 4      | 0      | 0     | 587   | 583   | 0.16  | 5.19     | 5.52     | 0.43
       0  | 235   | 646    | 0      | 0      | 0     | 563   | 560   | 0.17  | 6.04     | 6.40     | 0.47
       0  | 236   | 718    | 7      | 0      | 0     | 622   | 612   | 0.17  | 5.50     | 5.88     | 0.45
       0  | 237   | 695    | 6      | 0      | 0     | 613   | 607   | 0.17  | 5.55     | 5.86     | 0.46
       0  | 238   | 721    | 8      | 0      | 0     | 642   | 636   | 0.17  | 5.72     | 6.08     | 0.46
       0  | 239   | 691    | 14     | 0      | 0     | 619   | 610   | 0.16  | 5.25     | 5.65     | 0.45
       0  | 240   | 710    | 3      | 0      | 0     | 636   | 630   | 0.17  | 5.45     | 5.79     | 0.45
       0  | 241   | 740    | 7      | 0      | 0     | 670   | 663   | 0.16  | 5.41     | 5.72     | 0.44
       0  | 242   | 666    | 3      | 0      | 0     | 597   | 589   | 0.16  | 5.46     | 5.84     | 0.46
       0  | 243   | 685    | 10     | 0      | 0     | 614   | 609   | 0.17  | 6.11     | 6.43     | 0.46
       0  | 244   | 641    | 4      | 0      | 0     | 567   | 562   | 0.17  | 5.47     | 5.84     | 0.46
       0  | 245   | 688    | 49     | 0      | 0     | 628   | 617   | 0.17  | 5.67     | 6.13     | 0.45
       0  | 246   | 675    | 32     | 0      | 0     | 610   | 601   | 0.17  | 5.67     | 6.00     | 0.46
       0  | 247   | 673    | 12     | 0      | 0     | 596   | 589   | 0.17  | 6.10     | 6.48     | 0.47
       0  | 248   | 709    | 9      | 0      | 0     | 642   | 631   | 0.16  | 5.43     | 5.83     | 0.46
       0  | 249   | 685    | 4      | 0      | 0     | 587   | 581   | 0.17  | 5.85     | 6.20     | 0.45
       0  | 250   | 692    | 6      | 0      | 0     | 611   | 603   | 0.17  | 5.58     | 5.94     | 0.45
       0  | 251   | 717    | 6      | 0      | 0     | 635   | 631   | 0.16  | 5.92     | 6.26     | 0.44
       0  | 252   | 752    | 4      | 0      | 0     | 681   | 676   | 0.16  | 4.96     | 5.28     | 0.43
       0  | 253   | 692    | 0      | 0      | 0     | 618   | 614   | 0.17  | 5.34     | 5.70     | 0.47
       0  | 254   | 658    | 4      | 0      | 0     | 588   | 580   | 0.17  | 5.96     | 6.35     | 0.48
       0  | 255   | 689    | 8      | 0      | 0     | 614   | 605   | 0.17  | 5.93     | 6.30     | 0.44
       0  | 256   | 658    | 14     | 0      | 0     | 597   | 586   | 0.18  | 6.26     | 6.63     | 0.46
       0  | 257   | 636    | 28     | 0      | 0     | 579   | 568   | 0.16  | 5.01     | 5.48     | 0.44
       0  | 258   | 675    | 20     | 0      | 0     | 589   | 585   | 0.17  | 5.74     | 6.05     | 0.46
       0  | 259   | 666    | 16     | 0      | 0     | 578   | 573   | 0.16  | 5.45     | 5.83     | 0.44
       0  | 260   | 705    | 3      | 0      | 0     | 624   | 621   | 0.17  | 5.76     | 6.07     | 0.47
       0  | 261   | 722    | 19     | 0      | 0     | 646   | 638   | 0.17  | 5.46     | 5.81     | 0.45
       0  | 262   | 735    | 12     | 0      | 0     | 672   | 658   | 0.17  | 5.41     | 5.80     | 0.46
       0  | 263   | 713    | 11     | 0      | 0     | 626   | 621   | 0.16  | 5.52     | 5.89     | 0.45
       0  | 264   | 713    | 15     | 0      | 0     | 636   | 628   | 0.16  | 5.28     | 5.59     | 0.45
       0  | 265   | 646    | 21     | 0      | 0     | 579   | 570   | 0.17  | 5.93     | 6.31     | 0.44
       0  | 266   | 654    | 28     | 0      | 0     | 598   | 590   | 0.17  | 5.37     | 5.81     | 0.47
       0  | 267   | 668    | 17     | 0      | 0     | 603   | 595   | 0.18  | 5.97     | 6.37     | 0.47
       0  | 268   | 678    | 22     | 0      | 0     | 602   | 592   | 0.17  | 5.81     | 6.21     | 0.47
       0  | 269   | 665    | 11     | 0      | 0     | 590   | 586   | 0.17  | 5.39     | 5.74     | 0.45
       0  | 270   | 677    | 13     | 0      | 0     | 610   | 603   | 0.17  | 5.52     | 5.90     | 0.45
       0  | 271   | 695    | 0      | 0      | 0     | 605   | 602   | 0.16  | 4.77     | 5.12     | 0.42
       0  | 272   | 676    | 4      | 0      | 0     | 609   | 606   | 0.17  | 5.22     | 5.59     | 0.46
       0  | 273   | 688    | 5      | 0      | 0     | 626   | 622   | 0.17  | 5.61     | 5.98     | 0.47
       0  | 274   | 757    | 7      | 0      | 0     | 668   | 659   | 0.17  | 5.66     | 6.06     | 0.46
       0  | 275   | 750    | 10     | 0      | 0     | 667   | 659   | 0.18  | 6.01     | 6.45     | 0.48
       0  | 276   | 721    | 10     | 0      | 0     | 634   | 629   | 0.17  | 5.43     | 5.82     | 0.47
       0  | 277   | 691    | 3      | 0      | 0     | 597   | 590   | 0.17  | 5.91     | 6.33     | 0.49
       0  | 278   | 616    | 15     | 0      | 0     | 546   | 538   | 0.17  | 5.68     | 6.05     | 0.45
       0  | 279   | 658    | 12     | 0      | 0     | 593   | 582   | 0.17  | 5.66     | 6.12     | 0.46
       0  | 280   | 663    | 3      | 0      | 0     | 576   | 570   | 0.18  | 6.39     | 6.82     | 0.47
       0  | 281   | 699    | 22     | 0      | 0     | 635   | 625   | 0.17  | 5.23     | 5.64     | 0.47
       0  | 282   | 689    | 29     | 0      | 0     | 617   | 612   | 0.16  | 5.48     | 5.86     | 0.45
       0  | 283   | 720    | 21     | 0      | 0     | 630   | 624   | 0.18  | 6.02     | 6.41     | 0.48
       0  | 284   | 694    | 0      | 0      | 0     | 605   | 600   | 0.17  | 5.71     | 6.09     | 0.46
       0  | 285   | 693    | 9      | 0      | 0     | 619   | 612   | 0.18  | 6.04     | 6.43     | 0.46
       0  | 286   | 697    | 14     | 0      | 0     | 623   | 622   | 0.17  | 5.64     | 5.97     | 0.46
       0  | 287   | 718    | 3      | 0      | 0     | 648   | 644   | 0.17  | 5.44     | 5.81     | 0.46
       0  | 288   | 721    | 14     | 0      | 0     | 648   | 638   | 0.16  | 5.19     | 5.61     | 0.45
       0  | 289   | 709    | 7      | 0      | 0     | 637   | 629   | 0.17  | 5.55     | 5.97     | 0.47
       0  | 290   | 671    | 14     | 0      | 0     | 585   | 580   | 0.17  | 6.24     | 6.60     | 0.47
       0  | 291   | 634    | 14     | 0      | 0     | 563   | 559   | 0.18  | 5.82     | 6.20     | 0.47
       0  | 292   | 706    | 12     | 0      | 0     | 629   | 622   | 0.17  | 5.89     | 6.28     | 0.46
       0  | 293   | 651    | 16     | 0      | 0     | 574   | 571   | 0.27  | 5.65     | 6.09     | 0.47
       0  | 294   | 664    | 6      | 0      | 0     | 579   | 572   | 0.17  | 5.37     | 5.74     | 0.47
       0  | 295   | 704    | 26     | 0      | 0     | 623   | 615   | 0.18  | 5.75     | 6.20     | 0.49
       0  | 296   | 665    | 15     | 0      | 0     | 586   | 578   | 0.17  | 5.39     | 5.81     | 0.46
       0  | 297   | 647    | 22     | 0      | 0     | 565   | 558   | 0.17  | 5.24     | 5.64     | 0.45
       0  | 298   | 697    | 18     | 0      | 0     | 616   | 609   | 0.17  | 5.20     | 5.57     | 0.43
       0  | 299   | 732    | 11     | 0      | 0     | 651   | 645   | 0.17  | 6.34     | 6.72     | 0.48
       0  | 300   | 720    | 19     | 0      | 0     | 635   | 625   | 0.18  | 6.17     | 6.53     | 0.48
       0  | 301   | 730    | 6      | 0      | 0     | 655   | 650   | 0.17  | 5.38     | 5.75     | 0.45
       0  | 302   | 706    | 6      | 0      | 0     | 623   | 617   | 0.18  | 5.69     | 6.07     | 0.47
       0  | 303   | 707    | 12     | 0      | 0     | 636   | 631   | 0.16  | 5.75     | 6.12     | 0.45
       0  | 304   | 664    | 0      | 0      | 0     | 584   | 580   | 0.18  | 6.10     | 6.48     | 0.48
       0  | 305   | 677    | 10     | 0      | 0     | 618   | 616   | 0.18  | 5.94     | 6.30     | 0.46
       0  | 306   | 673    | 0      | 0      | 0     | 576   | 574   | 0.17  | 5.60     | 5.94     | 0.46
       0  | 307   | 674    | 3      | 0      | 0     | 589   | 584   | 0.19  | 6.34     | 6.70     | 0.48
       0  | 308   | 682    | 8      | 0      | 0     | 610   | 608   | 0.17  | 5.59     | 5.94     | 0.45
       0  | 309   | 664    | 7      | 0      | 0     | 575   | 573   | 0.17  | 5.38     | 5.79     | 0.47
       0  | 310   | 681    | 22     | 0      | 0     | 609   | 598   | 0.17  | 5.59     | 6.06     | 0.47
       0  | 311   | 694    | 14     | 0      | 0     | 641   | 630   | 0.16  | 5.15     | 5.58     | 0.45
       0  | 312   | 685    | 10     | 0      | 0     | 608   | 604   | 0.18  | 6.07     | 6.45     | 0.47
       0  | 313   | 737    | 26     | 0      | 0     | 656   | 651   | 0.17  | 5.85     | 6.23     | 0.47
       0  | 314   | 698    | 18     | 0      | 0     | 641   | 628   | 0.16  | 5.00     | 5.39     | 0.44
       0  | 315   | 703    | 21     | 0      | 0     | 625   | 622   | 0.17  | 5.26     | 5.61     | 0.44
       0  | 316   | 662    | 18     | 0      | 0     | 590   | 581   | 0.17  | 5.41     | 5.82     | 0.46
       0  | 317   | 709    | 17     | 0      | 0     | 629   | 625   | 0.16  | 5.14     | 5.48     | 0.44
       0  | 318   | 699    | 18     | 0      | 0     | 633   | 623   | 0.18  | 6.01     | 6.39     | 0.45
       0  | 319   | 654    | 10     | 0      | 0     | 575   | 566   | 0.18  | 5.87     | 6.30     | 0.48
       0  | 320   | 678    | 14     | 0      | 0     | 607   | 601   | 0.18  | 6.16     | 6.56     | 0.48
       0  | 321   | 687    | 14     | 0      | 0     | 606   | 598   | 0.18  | 6.08     | 6.45     | 0.47
       0  | 322   | 676    | 4      | 0      | 0     | 594   | 592   | 0.18  | 5.87     | 6.26     | 0.48
       0  | 323   | 672    | 0      | 0      | 0     | 591   | 586   | 0.18  | 5.46     | 5.84     | 0.44
       0  | 324   | 712    | 6      | 0      | 0     | 630   | 623   | 0.17  | 5.62     | 6.01     | 0.47
       0  | 325   | 692    | 8      | 0      | 0     | 604   | 599   | 0.18  | 6.20     | 6.57     | 0.46
       0  | 326   | 688    | 17     | 0      | 0     | 610   | 602   | 0.17  | 6.03     | 6.42     | 0.47
       0  | 327   | 687    | 17     | 0      | 0     | 618   | 610   | 0.17  | 5.91     | 6.32     | 0.47
       0  | 328   | 712    | 3      | 0      | 0     | 622   | 620   | 0.17  | 5.43     | 5.80     | 0.47
       0  | 329   | 721    | 31     | 0      | 0     | 653   | 644   | 0.17  | 5.32     | 5.74     | 0.45
       0  | 330   | 703    | 12     | 0      | 0     | 632   | 625   | 0.17  | 4.81     | 5.21     | 0.45
       0  | 331   | 680    | 20     | 0      | 0     | 625   | 615   | 0.17  | 5.15     | 5.61     | 0.46
       0  | 332   | 690    | 18     | 0      | 0     | 617   | 605   | 0.17  | 5.37     | 5.80     | 0.45
       0  | 333   | 679    | 9      | 0      | 0     | 610   | 605   | 0.17  | 5.64     | 6.01     | 0.47
       0  | 334   | 663    | 9      | 0      | 0     | 591   | 587   | 0.19  | 6.50     | 6.89     | 0.48
       0  | 335   | 696    | 9      | 0      | 0     | 612   | 604   | 0.17  | 5.36     | 5.77     | 0.44
       0  | 336   | 678    | 7      | 0      | 0     | 590   | 586   | 0.18  | 6.14     | 6.54     | 0.49
       0  | 337   | 659    | 3      | 0      | 0     | 571   | 564   | 0.18  | 5.58     | 5.95     | 0.46
       0  | 338   | 695    | 0      | 0      | 0     | 618   | 614   | 0.19  | 6.18     | 6.61     | 0.47
       0  | 339   | 709    | 3      | 0      | 0     | 619   | 609   | 0.17  | 5.41     | 5.81     | 0.46
       0  | 340   | 710    | 8      | 0      | 0     | 621   | 615   | 0.17  | 5.64     | 6.07     | 0.47
       0  | 341   | 681    | 14     | 0      | 0     | 593   | 584   | 0.17  | 5.56     | 5.98     | 0.46
       0  | 342   | 661    | 14     | 0      | 0     | 575   | 569   | 0.17  | 5.31     | 5.73     | 0.46
       0  | 343   | 721    | 6      | 0      | 0     | 643   | 639   | 0.17  | 5.17     | 5.58     | 0.46
       0  | 344   | 663    | 10     | 0      | 0     | 594   | 585   | 0.17  | 5.68     | 6.12     | 0.48
       0  | 345   | 698    | 6      | 0      | 0     | 617   | 610   | 0.17  | 5.57     | 5.99     | 0.46
       0  | 346   | 696    | 10     | 0      | 0     | 626   | 624   | 0.17  | 5.36     | 5.72     | 0.44
       0  | 347   | 678    | 12     | 0      | 0     | 605   | 599   | 0.17  | 5.66     | 6.07     | 0.47
       0  | 348   | 713    | 24     | 0      | 0     | 640   | 628   | 0.17  | 5.63     | 6.08     | 0.44
       0  | 349   | 725    | 14     | 0      | 0     | 648   | 644   | 0.19  | 6.32     | 6.68     | 0.49
       0  | 350   | 717    | 17     | 0      | 0     | 637   | 628   | 0.18  | 6.04     | 6.48     | 0.48
       0  | 351   | 679    | 13     | 0      | 0     | 598   | 589   | 0.17  | 5.43     | 5.81     | 0.45
       0  | 352   | 667    | 15     | 0      | 0     | 587   | 582   | 0.18  | 6.03     | 6.46     | 0.49
       0  | 353   | 693    | 22     | 0      | 0     | 616   | 602   | 0.18  | 5.96     | 6.40     | 0.47
       0  | 354   | 673    | 3      | 0      | 0     | 584   | 580   | 0.17  | 5.52     | 5.89     | 0.46
       0  | 355   | 670    | 3      | 0      | 0     | 608   | 605   | 0.17  | 5.54     | 5.95     | 0.46
       0  | 356   | 669    | 10     | 0      | 0     | 600   | 592   | 0.18  | 6.01     | 6.43     | 0.45
       0  | 357   | 682    | 18     | 0      | 0     | 603   | 597   | 0.18  | 5.76     | 6.13     | 0.45
       0  | 358   | 672    | 3      | 0      | 0     | 602   | 595   | 0.17  | 5.22     | 5.64     | 0.45
       0  | 359   | 680    | 16     | 0      | 0     | 619   | 612   | 0.17  | 5.89     | 6.30     | 0.47
       0  | 360   | 694    | 17     | 0      | 0     | 619   | 614   | 0.17  | 5.64     | 6.05     | 0.46
       0  | 361   | 688    | 3      | 0      | 0     | 603   | 600   | 0.19  | 6.73     | 7.09     | 0.49
       0  | 362   | 662    | 0      | 0      | 0     | 581   | 576   | 0.17  | 5.10     | 5.54     | 0.47
       0  | 363   | 726    | 3      | 0      | 0     | 642   | 634   | 0.18  | 5.66     | 6.08     | 0.48
       0  | 364   | 690    | 20     | 0      | 0     | 630   | 618   | 0.19  | 5.77     | 6.22     | 0.46
       0  | 365   | 724    | 40     | 0      | 0     | 658   | 645   | 0.17  | 5.10     | 5.55     | 0.46
       0  | 366   | 707    | 20     | 0      | 0     | 631   | 624   | 0.18  | 6.09     | 6.47     | 0.47
       0  | 367   | 698    | 11     | 0      | 0     | 631   | 622   | 0.18  | 5.21     | 5.66     | 0.46
       0  | 368   | 716    | 6      | 0      | 0     | 636   | 632   | 0.18  | 5.82     | 6.25     | 0.48
       0  | 369   | 711    | 17     | 0      | 0     | 628   | 621   | 0.18  | 5.76     | 6.19     | 0.49
       0  | 370   | 692    | 36     | 0      | 0     | 623   | 615   | 0.17  | 5.78     | 6.20     | 0.44
       0  | 371   | 656    | 10     | 0      | 0     | 586   | 581   | 0.17  | 5.84     | 6.25     | 0.47
       0  | 372   | 638    | 7      | 0      | 0     | 560   | 559   | 0.18  | 6.84     | 7.16     | 0.47
       0  | 373   | 672    | 0      | 0      | 0     | 596   | 590   | 0.17  | 5.77     | 6.16     | 0.46
       0  | 374   | 706    | 0      | 0      | 0     | 628   | 628   | 0.18  | 5.83     | 6.21     | 0.48
       0  | 375   | 700    | 10     | 0      | 0     | 595   | 590   | 0.17  | 5.22     | 5.65     | 0.45
       0  | 376   | 686    | 10     | 0      | 0     | 607   | 600   | 0.18  | 5.33     | 5.76     | 0.46
       0  | 377   | 680    | 11     | 0      | 0     | 612   | 607   | 0.19  | 6.21     | 6.63     | 0.48
       0  | 378   | 708    | 7      | 0      | 0     | 617   | 616   | 0.17  | 5.37     | 5.75     | 0.47
       0  | 379   | 678    | 18     | 0      | 0     | 597   | 586   | 0.17  | 5.88     | 6.31     | 0.47
       0  | 380   | 687    | 16     | 0      | 0     | 608   | 602   | 0.18  | 5.47     | 5.90     | 0.46
       0  | 381   | 656    | 27     | 0      | 0     | 595   | 590   | 0.17  | 4.86     | 5.29     | 0.45
       0  | 382   | 705    | 12     | 0      | 0     | 624   | 618   | 0.18  | 5.81     | 6.24     | 0.47
       0  | 383   | 628    | 17     | 0      | 0     | 564   | 553   | 0.17  | 5.09     | 5.55     | 0.46
       0  | 384   | 691    | 15     | 0      | 0     | 626   | 614   | 0.17  | 5.28     | 5.77     | 0.46
       0  | 385   | 748    | 3      | 0      | 0     | 665   | 658   | 0.18  | 5.26     | 5.72     | 0.49
       0  | 386   | 708    | 9      | 0      | 0     | 629   | 624   | 0.18  | 5.40     | 5.76     | 0.46
       0  | 387   | 658    | 6      | 0      | 0     | 580   | 575   | 0.19  | 6.12     | 6.53     | 0.48
       0  | 388   | 745    | 3      | 0      | 0     | 683   | 677   | 0.18  | 5.69     | 6.08     | 0.48
       0  | 389   | 682    | 3      | 0      | 0     | 602   | 595   | 0.19  | 6.16     | 6.59     | 0.48
       0  | 390   | 725    | 8      | 0      | 0     | 637   | 635   | 0.18  | 5.90     | 6.29     | 0.48
       0  | 391   | 688    | 8      | 0      | 0     | 613   | 610   | 0.19  | 6.71     | 7.10     | 0.48
       0  | 392   | 703    | 4      | 0      | 0     | 617   | 614   | 0.18  | 6.16     | 6.53     | 0.47
       0  | 393   | 697    | 21     | 0      | 0     | 620   | 613   | 0.18  | 5.94     | 6.35     | 0.46
       0  | 394   | 643    | 27     | 0      | 0     | 586   | 574   | 0.18  | 5.45     | 5.94     | 0.45
       0  | 395   | 679    | 28     | 0      | 0     | 590   | 582   | 0.18  | 5.86     | 6.28     | 0.46
       0  | 396   | 700    | 13     | 0      | 0     | 601   | 594   | 0.18  | 5.65     | 6.10     | 0.48
       0  | 397   | 672    | 20     | 0      | 0     | 613   | 602   | 0.18  | 5.58     | 6.05     | 0.48
       0  | 398   | 700    | 15     | 0      | 0     | 626   | 614   | 0.18  | 5.21     | 5.67     | 0.47
       0  | 399   | 680    | 23     | 0      | 0     | 621   | 611   | 0.18  | 5.29     | 5.81     | 0.47
       0  | 400   | 670    | 22     | 0      | 0     | 614   | 606   | 0.18  | 5.46     | 5.86     | 0.45
       0  | 401   | 678    | 18     | 0      | 0     | 615   | 605   | 0.17  | 5.21     | 5.64     | 0.45
       0  | 402   | 662    | 3      | 0      | 0     | 573   | 569   | 0.17  | 5.16     | 5.57     | 0.46
       0  | 403   | 662    | 6      | 0      | 0     | 607   | 604   | 0.18  | 6.11     | 6.50     | 0.46
       0  | 404   | 637    | 8      | 0      | 0     | 544   | 537   | 0.17  | 5.09     | 5.51     | 0.45
       0  | 405   | 756    | 17     | 0      | 0     | 675   | 669   | 0.18  | 6.33     | 6.75     | 0.46
       0  | 406   | 696    | 4      | 0      | 0     | 612   | 605   | 0.18  | 5.74     | 6.23     | 0.51
       0  | 407   | 726    | 4      | 0      | 0     | 640   | 637   | 0.19  | 5.76     | 6.21     | 0.49
       0  | 408   | 632    | 6      | 0      | 0     | 563   | 562   | 0.18  | 5.83     | 6.18     | 0.46
       0  | 409   | 735    | 7      | 0      | 0     | 645   | 639   | 0.18  | 5.19     | 5.68     | 0.47
       0  | 410   | 739    | 10     | 0      | 0     | 659   | 653   | 0.18  | 5.02     | 5.47     | 0.44
       0  | 411   | 641    | 3      | 0      | 0     | 566   | 563   | 0.19  | 6.54     | 6.95     | 0.50
       0  | 412   | 674    | 8      | 0      | 0     | 587   | 577   | 0.18  | 6.21     | 6.63     | 0.46
       0  | 413   | 713    | 18     | 0      | 0     | 636   | 631   | 0.17  | 5.28     | 5.73     | 0.45
       0  | 414   | 687    | 25     | 0      | 0     | 619   | 609   | 0.18  | 5.47     | 5.91     | 0.46
       0  | 415   | 710    | 21     | 0      | 0     | 643   | 638   | 0.19  | 5.26     | 5.68     | 0.47
       0  | 416   | 699    | 22     | 0      | 0     | 627   | 614   | 0.18  | 5.52     | 5.99     | 0.48
       0  | 417   | 690    | 17     | 0      | 0     | 616   | 610   | 0.18  | 5.59     | 6.06     | 0.47
       0  | 418   | 682    | 11     | 0      | 0     | 596   | 592   | 0.18  | 5.43     | 5.88     | 0.46
       0  | 419   | 720    | 14     | 0      | 0     | 635   | 628   | 0.19  | 5.92     | 6.34     | 0.46
       0  | 420   | 725    | 12     | 0      | 0     | 640   | 634   | 0.19  | 6.36     | 6.78     | 0.47
       0  | 421   | 689    | 0      | 0      | 0     | 601   | 600   | 0.19  | 5.61     | 6.03     | 0.49
       0  | 422   | 675    | 9      | 0      | 0     | 617   | 612   | 0.19  | 6.07     | 6.47     | 0.47
       0  | 423   | 651    | 9      | 0      | 0     | 569   | 565   | 0.18  | 5.18     | 5.56     | 0.44
       0  | 424   | 691    | 0      | 0      | 0     | 596   | 594   | 0.17  | 4.86     | 5.30     | 0.46
       0  | 425   | 708    | 0      | 0      | 0     | 627   | 625   | 0.18  | 4.83     | 5.25     | 0.46
       0  | 426   | 692    | 9      | 0      | 0     | 600   | 592   | 0.18  | 5.30     | 5.72     | 0.47
       0  | 427   | 700    | 7      | 0      | 0     | 624   | 614   | 0.18  | 5.42     | 5.89     | 0.46
       0  | 428   | 650    | 7      | 0      | 0     | 582   | 578   | 0.18  | 4.74     | 5.16     | 0.44
       0  | 429   | 657    | 20     | 0      | 0     | 588   | 579   | 0.17  | 5.15     | 5.59     | 0.45
       0  | 430   | 715    | 12     | 0      | 0     | 637   | 629   | 0.18  | 5.29     | 5.70     | 0.46
       0  | 431   | 636    | 15     | 0      | 0     | 578   | 571   | 0.18  | 6.19     | 6.64     | 0.48
       0  | 432   | 686    | 12     | 0      | 0     | 625   | 616   | 0.18  | 5.21     | 5.70     | 0.47
       0  | 433   | 699    | 21     | 0      | 0     | 640   | 628   | 0.18  | 5.90     | 6.38     | 0.46
       0  | 434   | 712    | 3      | 0      | 0     | 639   | 638   | 0.18  | 5.19     | 5.56     | 0.44
       0  | 435   | 695    | 21     | 0      | 0     | 622   | 616   | 0.18  | 5.40     | 5.82     | 0.46
       0  | 436   | 669    | 20     | 0      | 0     | 595   | 589   | 0.17  | 5.02     | 5.48     | 0.46
       0  | 437   | 682    | 21     | 0      | 0     | 609   | 599   | 0.18  | 5.27     | 5.71     | 0.46
       0  | 438   | 658    | 10     | 0      | 0     | 589   | 583   | 0.18  | 5.51     | 5.92     | 0.46
       0  | 439   | 672    | 16     | 0      | 0     | 606   | 600   | 0.18  | 5.66     | 6.09     | 0.48
       0  | 440   | 684    | 0      | 0      | 0     | 581   | 573   | 0.18  | 5.38     | 5.82     | 0.46
       0  | 441   | 734    | 3      | 0      | 0     | 647   | 644   | 0.18  | 5.54     | 6.00     | 0.47
       0  | 442   | 719    | 4      | 0      | 0     | 621   | 612   | 0.18  | 5.78     | 6.20     | 0.45
       0  | 443   | 683    | 5      | 0      | 0     | 591   | 588   | 0.18  | 5.00     | 5.47     | 0.46
       0  | 444   | 687    | 10     | 0      | 0     | 608   | 600   | 0.19  | 5.83     | 6.30     | 0.48
       0  | 445   | 689    | 22     | 0      | 0     | 624   | 618   | 0.18  | 5.65     | 6.04     | 0.46
       0  | 446   | 673    | 13     | 0      | 0     | 582   | 577   | 0.18  | 5.24     | 5.73     | 0.48
       0  | 447   | 685    | 17     | 0      | 0     | 616   | 608   | 0.18  | 4.82     | 5.27     | 0.45
       0  | 448   | 679    | 12     | 0      | 0     | 619   | 612   | 0.18  | 5.74     | 6.13     | 0.45
       0  | 449   | 656    | 6      | 0      | 0     | 589   | 582   | 0.18  | 5.15     | 5.60     | 0.45
       0  | 450   | 699    | 14     | 0      | 0     | 618   | 610   | 0.18  | 5.47     | 5.95     | 0.47
       0  | 451   | 737    | 6      | 0      | 0     | 651   | 642   | 0.19  | 6.05     | 6.51     | 0.47
       0  | 452   | 736    | 11     | 0      | 0     | 652   | 648   | 0.18  | 5.51     | 5.89     | 0.45
       0  | 453   | 690    | 6      | 0      | 0     | 610   | 608   | 0.18  | 4.79     | 5.24     | 0.46
       0  | 454   | 640    | 8      | 0      | 0     | 571   | 567   | 0.19  | 5.26     | 5.78     | 0.48
       0  | 455   | 668    | 19     | 0      | 0     | 613   | 608   | 0.19  | 5.70     | 6.17     | 0.47
       0  | 456   | 659    | 11     | 0      | 0     | 576   | 571   | 0.18  | 5.53     | 5.97     | 0.46
       0  | 457   | 739    | 4      | 0      | 0     | 660   | 654   | 0.18  | 5.34     | 5.83     | 0.47
       0  | 458   | 711    | 34     | 0      | 0     | 619   | 612   | 0.17  | 4.94     | 5.46     | 0.46
       0  | 459   | 667    | 31     | 0      | 0     | 579   | 573   | 0.18  | 5.66     | 6.14     | 0.48
       0  | 460   | 691    | 12     | 0      | 0     | 605   | 598   | 0.19  | 5.86     | 6.36     | 0.48
       0  | 461   | 685    | 13     | 0      | 0     | 613   | 607   | 0.19  | 6.57     | 6.99     | 0.47
       0  | 462   | 718    | 6      | 0      | 0     | 643   | 631   | 0.19  | 5.95     | 6.44     | 0.48
       0  | 463   | 694    | 34     | 0      | 0     | 620   | 611   | 0.19  | 6.20     | 6.68     | 0.48
       0  | 464   | 703    | 11     | 0      | 0     | 616   | 612   | 0.18  | 5.54     | 6.00     | 0.47
       0  | 465   | 675    | 8      | 0      | 0     | 606   | 603   | 0.19  | 5.02     | 5.43     | 0.47
       0  | 466   | 684    | 9      | 0      | 0     | 603   | 589   | 0.19  | 5.83     | 6.32     | 0.46
       0  | 467   | 691    | 24     | 0      | 0     | 621   | 612   | 0.18  | 5.51     | 5.98     | 0.48
       0  | 468   | 699    | 19     | 0      | 0     | 632   | 627   | 0.19  | 5.45     | 5.96     | 0.48
       0  | 469   | 675    | 15     | 0      | 0     | 614   | 605   | 0.19  | 5.62     | 6.12     | 0.46
       0  | 470   | 723    | 30     | 0      | 0     | 654   | 644   | 0.18  | 5.08     | 5.61     | 0.46
       0  | 471   | 648    | 4      | 0      | 0     | 571   | 563   | 0.18  | 5.10     | 5.62     | 0.47
       0  | 472   | 691    | 16     | 0      | 0     | 624   | 619   | 0.32  | 5.59     | 6.18     | 0.47
       0  | 473   | 731    | 35     | 0      | 0     | 627   | 621   | 0.18  | 5.12     | 5.57     | 0.46
       0  | 474   | 652    | 6      | 0      | 0     | 594   | 591   | 0.19  | 5.91     | 6.36     | 0.48
       0  | 475   | 702    | 12     | 0      | 0     | 615   | 610   | 0.18  | 5.08     | 5.53     | 0.46
       0  | 476   | 729    | 11     | 0      | 0     | 636   | 624   | 0.19  | 5.24     | 5.81     | 0.49
       0  | 477   | 677    | 13     | 0      | 0     | 593   | 587   | 0.20  | 6.07     | 6.54     | 0.49
       0  | 478   | 677    | 6      | 0      | 0     | 583   | 580   | 0.19  | 5.41     | 5.88     | 0.46
       0  | 479   | 686    | 20     | 0      | 0     | 612   | 604   | 0.18  | 4.66     | 5.15     | 0.44
       0  | 480   | 703    | 16     | 0      | 0     | 618   | 614   | 0.20  | 5.78     | 6.23     | 0.48
       0  | 481   | 682    | 6      | 0      | 0     | 608   | 598   | 0.18  | 5.44     | 5.91     | 0.46
       0  | 482   | 677    | 18     | 0      | 0     | 600   | 596   | 0.18  | 4.97     | 5.43     | 0.45
       0  | 483   | 756    | 11     | 0      | 0     | 673   | 669   | 0.19  | 5.98     | 6.43     | 0.47
       0  | 484   | 627    | 11     | 0      | 0     | 563   | 555   | 0.19  | 5.21     | 5.73     | 0.47
       0  | 485   | 680    | 34     | 0      | 0     | 618   | 610   | 0.20  | 5.75     | 6.25     | 0.48
       0  | 486   | 699    | 15     | 0      | 0     | 636   | 629   | 0.18  | 4.86     | 5.34     | 0.46
       0  | 487   | 679    | 17     | 0      | 0     | 610   | 600   | 0.20  | 5.77     | 6.27     | 0.49
       0  | 488   | 690    | 45     | 0      | 0     | 622   | 615   | 0.18  | 5.22     | 5.67     | 0.45
       0  | 489   | 724    | 4      | 0      | 0     | 637   | 631   | 0.18  | 5.34     | 5.80     | 0.45
       0  | 490   | 686    | 19     | 0      | 0     | 630   | 618   | 0.19  | 5.39     | 5.96     | 0.47
       0  | 491   | 650    | 14     | 0      | 0     | 571   | 566   | 0.19  | 5.50     | 5.96     | 0.48
       0  | 492   | 702    | 14     | 0      | 0     | 615   | 604   | 0.19  | 5.73     | 6.22     | 0.45
       0  | 493   | 684    | 7      | 0      | 0     | 599   | 591   | 0.19  | 5.73     | 6.23     | 0.47
       0  | 494   | 685    | 23     | 0      | 0     | 605   | 599   | 0.20  | 5.92     | 6.41     | 0.49
       0  | 495   | 684    | 11     | 0      | 0     | 622   | 610   | 0.18  | 5.06     | 5.54     | 0.45
       0  | 496   | 721    | 23     | 0      | 0     | 646   | 640   | 0.19  | 5.76     | 6.23     | 0.46
       0  | 497   | 691    | 7      | 0      | 0     | 604   | 601   | 0.19  | 5.61     | 6.03     | 0.47
       0  | 498   | 656    | 37     | 0      | 0     | 591   | 581   | 0.19  | 5.62     | 6.03     | 0.46
       0  | 499   | 707    | 27     | 0      | 0     | 626   | 615   | 0.19  | 5.57     | 6.14     | 0.47
       0  | 500   | 671    | 15     | 0      | 0     | 612   | 601   | 0.19  | 5.07     | 5.62     | 0.48
       0  | 501   | 690    | 25     | 0      | 0     | 631   | 618   | 0.19  | 5.55     | 6.09     | 0.46
       0  | 502   | 666    | 6      | 0      | 0     | 576   | 575   | 0.18  | 5.06     | 5.56     | 0.47
       0  | 503   | 714    | 20     | 0      | 0     | 647   | 636   | 0.19  | 6.12     | 6.58     | 0.46
       0  | 504   | 711    | 14     | 0      | 0     | 641   | 638   | 0.20  | 6.36     | 6.72     | 0.45
       0  | 505   | 653    | 13     | 0      | 0     | 573   | 565   | 0.19  | 5.65     | 6.16     | 0.48
       0  | 506   | 698    | 8      | 0      | 0     | 635   | 631   | 0.19  | 5.63     | 6.10     | 0.46
       0  | 507   | 690    | 17     | 0      | 0     | 601   | 597   | 0.19  | 5.83     | 6.31     | 0.49
       0  | 508   | 695    | 3      | 0      | 0     | 623   | 618   | 0.19  | 5.53     | 6.07     | 0.48
       0  | 509   | 637    | 12     | 0      | 0     | 561   | 555   | 0.19  | 5.35     | 5.88     | 0.46
       0  | 510   | 706    | 4      | 0      | 0     | 617   | 610   | 0.18  | 5.19     | 5.60     | 0.44
       0  | 511   | 736    | 18     | 0      | 0     | 645   | 637   | 0.19  | 5.39     | 5.89     | 0.48
       0  | 512   | 690    | 20     | 0      | 0     | 602   | 593   | 0.19  | 5.21     | 5.69     | 0.46
       0  | 513   | 696    | 14     | 0      | 0     | 622   | 618   | 0.19  | 5.51     | 5.96     | 0.46
       0  | 514   | 713    | 16     | 0      | 0     | 649   | 641   | 0.19  | 5.90     | 6.36     | 0.46
       0  | 515   | 666    | 13     | 0      | 0     | 609   | 602   | 0.20  | 5.43     | 5.97     | 0.51
       0  | 516   | 642    | 22     | 0      | 0     | 583   | 572   | 0.19  | 5.89     | 6.38     | 0.49
       0  | 517   | 697    | 9      | 0      | 0     | 623   | 615   | 0.18  | 4.72     | 5.21     | 0.47
       0  | 518   | 723    | 15     | 0      | 0     | 639   | 635   | 0.19  | 5.15     | 5.63     | 0.47
       0  | 519   | 730    | 19     | 0      | 0     | 642   | 637   | 0.20  | 5.30     | 5.75     | 0.47
       0  | 520   | 679    | 3      | 0      | 0     | 612   | 610   | 0.19  | 5.22     | 5.64     | 0.45
       0  | 521   | 661    | 7      | 0      | 0     | 584   | 576   | 0.19  | 5.67     | 6.16     | 0.46
       0  | 522   | 676    | 8      | 0      | 0     | 589   | 585   | 0.19  | 5.33     | 5.83     | 0.47
       0  | 523   | 709    | 3      | 0      | 0     | 616   | 609   | 0.19  | 5.25     | 5.73     | 0.48
       0  | 524   | 660    | 9      | 0      | 0     | 569   | 565   | 0.19  | 5.52     | 6.00     | 0.47
       0  | 525   | 743    | 11     | 0      | 0     | 649   | 643   | 0.20  | 5.74     | 6.26     | 0.48
       0  | 526   | 680    | 1      | 0      | 0     | 590   | 588   | 0.19  | 5.39     | 5.87     | 0.47
       0  | 527   | 691    | 16     | 0      | 0     | 608   | 603   | 0.19  | 5.55     | 6.07     | 0.49
       0  | 528   | 682    | 8      | 0      | 0     | 607   | 606   | 0.19  | 5.21     | 5.69     | 0.46
       0  | 529   | 728    | 5      | 0      | 0     | 652   | 648   | 0.19  | 5.60     | 6.07     | 0.47
       0  | 530   | 684    | 18     | 0      | 0     | 601   | 596   | 0.20  | 5.56     | 6.02     | 0.47
       0  | 531   | 652    | 11     | 0      | 0     | 588   | 586   | 0.19  | 4.95     | 5.42     | 0.47
       0  | 532   | 691    | 9      | 0      | 0     | 635   | 632   | 0.20  | 6.18     | 6.66     | 0.49
       0  | 533   | 656    | 15     | 0      | 0     | 593   | 589   | 0.19  | 4.78     | 5.26     | 0.46
       0  | 534   | 712    | 16     | 0      | 0     | 641   | 636   | 0.19  | 5.47     | 5.94     | 0.48
       0  | 535   | 681    | 24     | 0      | 0     | 634   | 624   | 0.20  | 5.36     | 5.93     | 0.48
       0  | 536   | 701    | 44     | 0      | 0     | 636   | 634   | 0.19  | 4.96     | 5.44     | 0.48
       0  | 537   | 601    | 88     | 0      | 0     | 615   | 606   | 0.19  | 5.31     | 5.79     | 0.46
       0  | 538   | 294    | 391    | 0      | 0     | 595   | 591   | 0.19  | 5.62     | 5.98     | 0.44
       0  | 539   | 0      | 663    | 0      | 0     | 575   | 570   | 0.19  | 5.20     | 5.53     | 0.48
       ------------------------------------------------------------------------------------------------------

       Summary vs resolution
       ------------------------------------------------------------------------------------------------------
       ID | d min | # full | # part | # over | # ice | # sum | # prf | <Ibg> | <I/sigI> | <I/sigI> | <CC prf>
          |       |        |        |        |       |       |       |       |  (sum)   |  (prf)   |
       ------------------------------------------------------------------------------------------------------
       0  | 1.17  | 353    | 3      | 0      | 0     | 264   | 199   | 0.04  | 0.38     | 0.53     | 0.09
       0  | 1.19  | 1124   | 6      | 0      | 0     | 998   | 872   | 0.04  | 0.45     | 0.53     | 0.08
       0  | 1.21  | 2480   | 13     | 0      | 0     | 2216  | 2019  | 0.05  | 0.51     | 0.59     | 0.09
       0  | 1.23  | 4061   | 27     | 0      | 0     | 3658  | 3479  | 0.05  | 0.55     | 0.69     | 0.11
       0  | 1.26  | 5880   | 33     | 0      | 0     | 5272  | 5054  | 0.05  | 0.58     | 0.76     | 0.12
       0  | 1.28  | 8009   | 46     | 0      | 0     | 7104  | 6853  | 0.06  | 0.65     | 0.84     | 0.14
       0  | 1.31  | 10557  | 62     | 0      | 0     | 9365  | 9113  | 0.06  | 0.77     | 0.98     | 0.17
       0  | 1.35  | 13778  | 81     | 0      | 0     | 12256 | 11949 | 0.07  | 0.91     | 1.14     | 0.20
       0  | 1.38  | 18760  | 106    | 0      | 0     | 16735 | 16377 | 0.07  | 0.98     | 1.23     | 0.22
       0  | 1.42  | 22406  | 164    | 0      | 0     | 20142 | 19929 | 0.08  | 1.20     | 1.47     | 0.25
       0  | 1.47  | 26201  | 688    | 0      | 0     | 23701 | 23387 | 0.09  | 1.45     | 1.76     | 0.28
       0  | 1.52  | 27814  | 879    | 0      | 0     | 24522 | 24386 | 0.09  | 1.75     | 2.09     | 0.33
       0  | 1.58  | 28063  | 820    | 0      | 0     | 25790 | 25674 | 0.10  | 2.16     | 2.53     | 0.39
       0  | 1.65  | 28167  | 785    | 0      | 0     | 24642 | 24533 | 0.12  | 2.67     | 3.08     | 0.45
       0  | 1.74  | 28310  | 751    | 0      | 0     | 25129 | 25016 | 0.14  | 3.48     | 3.92     | 0.52
       0  | 1.85  | 28421  | 875    | 0      | 0     | 25992 | 25870 | 0.18  | 4.80     | 5.29     | 0.59
       0  | 1.99  | 28485  | 934    | 0      | 0     | 25023 | 24917 | 0.24  | 6.51     | 7.04     | 0.66
       0  | 2.19  | 28821  | 807    | 0      | 0     | 26107 | 26007 | 0.28  | 8.75     | 9.26     | 0.69
       0  | 2.51  | 28910  | 806    | 0      | 0     | 25603 | 25499 | 0.34  | 12.61    | 13.06    | 0.71
       0  | 3.17  | 29208  | 1087   | 0      | 0     | 26389 | 26267 | 0.41  | 25.23    | 25.29    | 0.72
       ------------------------------------------------------------------------------------------------------

       Summary for experiment 0
       ----------------------------------------------------------------
       Item                                  | Overall | Low    | High
       ----------------------------------------------------------------
       dmin                                  | 1.17    | 1.47   | 1.17
       dmax                                  | 151.27  | 151.27 | 1.47
       number fully recorded                 | 369808  | 282400 | 87408
       number partially recorded             | 8973    | 8432   | 541
       number with invalid background pixels | 70202   | 51241  | 18961
       number with invalid foreground pixels | 45552   | 35613  | 9939
       number with overloaded pixels         | 0       | 0      | 0
       number in powder rings                | 0       | 0      | 0
       number processed with summation       | 330908  | 252898 | 78010
       number processed with profile fitting | 327400  | 251556 | 75844
       number failed in background modelling | 7774    | 7321   | 453
       number failed in summation            | 45552   | 35613  | 9939
       number failed in profile fitting      | 49060   | 36955  | 12105
       <ibg>                                 | 0.17    | 0.20   | 0.07
       <i/sigi> (summation)                  | 5.62    | 7.08   | 0.90
       <i/sigi> (profile fitting)            | 6.01    | 7.48   | 1.13
       <cc prf>                              | 0.46    | 0.54   | 0.19
       cc_pearson sum/prf                    | 1.00    | 1.00   | 0.89
       cc_spearman sum/prf                   | 0.96    | 0.98   | 0.77
       ----------------------------------------------------------------

       ----------------------------------
               Read time | 140.77 seconds
            Extract time |   3.52 seconds
        Pre-process time |   0.47 seconds
            Process time | 226.68 seconds
       Post-process time |   0.00 seconds
              Total time | 376.31 seconds
               User time |   0.00 seconds
       ----------------------------------

      Saving 378781 reflections to integrated.pickle
       time taken: 0.522675
      Saving the experiments to integrated_experiments.json
       time taken: 0.0149889

      Total time taken: 161.655426


Checking this output we see that after loading in the reference reflections
from :file:`refined.pickle`,
new predictions are made up to the highest resolution at the corner of the
detector. This is fine, but if we wanted to we could have adjusted the
resolution limits using parameters :samp:`dmin` and :samp:`dmax`. The predictions
are made using the scan-varying crystal model recorded in
:file:`refined_experiments.json`. This ensures that prediction is made using
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

  hklout = "integrated.mtz"
  input {
    experiments = refined_experiments.json
    reflections = integrated.pickle
  }

  Removing 23949 reflections with negative variance
  Removing 27432 profile reflections with negative variance
  Removing 2 reflections with I/Sig(I) < -5.0
  Removing 0 profile reflections with I/Sig(I) < -5.0
  Removing 4039 incomplete reflections
  Title: from dials.export_mtz
  Space group symbol from file: P4
  Space group number from file: 75
  Space group from matrices: P 4 (No. 75)
  Point group symbol from file: 4
  Number of batches: 540
  Number of crystals: 1
  Number of Miller indices: 323359
  Resolution range: 150.012 1.17004
  History:
  Crystal 1:
    Name: XTAL
    Project: DIALS
    Id: 1
    Unit cell: (57.7873, 57.7873, 150.012, 90, 90, 90)
    Number of datasets: 1
    Dataset 1:
      Name: FROMDIALS
      Id: 1
      Wavelength: 0.97625
      Number of columns: 15
      label        #valid  %valid    min     max type
      H            323359 100.00%   0.00   46.00 H: index h,k,l
      K            323359 100.00%   0.00   47.00 H: index h,k,l
      L            323359 100.00%   0.00  114.00 H: index h,k,l
      M_ISYM       323359 100.00%   1.00    8.00 Y: M/ISYM, packed partial/reject flag and symmetry number
      BATCH        323359 100.00%   2.00  539.00 B: BATCH number
      IPR          323359 100.00%  -2.22 3937.64 J: intensity
      SIGIPR       323359 100.00%   0.05   62.78 Q: standard deviation
      I            323359 100.00% -32.84 4228.04 J: intensity
      SIGI         323359 100.00%   0.11   65.18 Q: standard deviation
      FRACTIONCALC 323359 100.00%   1.00    1.00 R: real
      XDET         323359 100.00%   6.60 2456.31 R: real
      YDET         323359 100.00%   5.79 2520.56 R: real
      ROT          323359 100.00%  82.01  162.69 R: real
      LP           323359 100.00%   0.00    0.76 R: real
      DQE          323359 100.00%   0.71    0.86 R: real



What to do Next
---------------

The following demonstrates how to take the output of dials processing and
continue with downstream analysis using the CCP4 programs pointless_, to sort the data and assign
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
  High resolution limit                       1.30      7.12      1.30

  Rmerge  (within I+/I-)                     0.061     0.024     0.416
  Rmerge  (all I+ and I-)                    0.069     0.026     0.488
  Rmeas (within I+/I-)                       0.075     0.030     0.575
  Rmeas (all I+ & I-)                        0.076     0.030     0.610
  Rpim (within I+/I-)                        0.043     0.017     0.395
  Rpim (all I+ & I-)                         0.033     0.014     0.358
  Rmerge in top intensity bin                0.029        -         -
  Total number of observations              308123      2257      5493
  Total number unique                        62352       499      2474
  Mean((I)/sd(I))                             10.7      27.1       1.4
  Mn(I) half-set correlation CC(1/2)         0.999     0.999     0.722
  Completeness                                98.2      99.8      80.1
  Multiplicity                                 4.9       4.5       2.2

  Anomalous completeness                      92.3     100.0      47.8
  Anomalous multiplicity                       2.4       3.0       1.5
  DelAnom correlation between half-sets     -0.002     0.279     0.065
  Mid-Slope of Anom Normal Probability       0.953       -         -


.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html
