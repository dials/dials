Advanced Tutorial
=================

Introduction
------------

DIALS processing may be performed by either running the individual tools (spot
finding, indexing, refinement, integration, exporting to MTZ) or you can run the
whole lot through :doc:`dials.process </programs/dials_process>`, which just
chains them together (and incidentally does all of the processing in P1.)

dials.process
-------------

In the simplest case, :doc:`dials.process </programs/dials_process>`
/here/are/all/images*.cbf` will do sensible processing, with a static model
of the experiment and sample, and will output a reflection file integrated.mtz
containing the intensity measurements assuming everything works correctly.
Some sensible options to use are:

 - :samp:`scan_varying=true` - allow the crystal orientation and unit cell
   constants to vary during the scan
 - :samp:`mp.nproc=1` - only use one processor (necessary currently for data in
   NeXus files)
 - :samp:`intensity.algorithm=sum` - use summation integtration, other
   algorithms are being added
 - :samp:`block_size=N` - for some N, split the data set into N degree blocks
   for integration, so as not to overload the computer
 - :samp:`-i` - pass the images to process through the standard input e.g. from
   :samp:`find . -name *.cbf` to avoid issues with limited command-line lengths

Running the Individual Steps: Macromolecule
-------------------------------------------

The following example uses a Thaumatin dataset collected using beamline I04
at Diamond Light Source which is available for download from |thaumatin|.

.. |thaumatin| image:: https://zenodo.org/badge/doi/10.5281/zenodo.10271.png
               :target: http://dx.doi.org/10.5281/zenodo.10271

A complete example script can be found
:download:`here<../user-tutorial/tutorial.sh>`, which can be run as follows::

  ./tutorial.sh /path/to/data

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that the image headers are read, and a file describing
their contents (:samp:`datablock.json`) is written. It's worth noting that if
this file is changed subsequent processing (even with :samp:`dials.process`) can
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
Here we tweak the minimum spot size (min_spot_size=3) and use multiple
processors to speed up the spot-finding (nproc=4).

::

  dials.find_spots datablock.json min_spot_size=3 nproc=4

This will just report the number of spots found.

::

  The following parameters have been modified:

  spotfinder {
    mp {
      nproc = 8
    }
    filter {
      min_spot_size = 3
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
  Extracted strong pixels from images.......................................21.09s
  Merged 8 pixel lists with 922120 pixels....................................0.02s
  Extracted 219125 spots.....................................................0.88s
  Calculated 219125 spot centroids...........................................0.79s
  Calculated 219125 spot intensities.........................................0.02s
  Filtered 116321 spots by number of pixels..................................0.01s
  Filtered 116082 spots by peak-centroid distance............................0.05s

  --------------------------------------------------------------------------------
  Saved 116082 reflections to strong.pickle..................................0.23s


Indexing
^^^^^^^^

The next step will be indexing of the strong spots, by default using a 3D FFT
algorithm, although the 1D FFT algorithm can be selected using the parameter
:samp:`indexing.method=fft1d`.

::

  dials.index datablock.json strong.pickle refinement.reflections.use_all_reflections=true

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

  Found max_cell: 229.7 Angstrom
  Setting d_min: 4.48575618871
  FFT gridding: (256,256,256)
  Number of centroids used: 8627
  model 1 (7863 reflections):
  Crystal:
      Unit cell: (58.179, 58.461, 149.622, 90.337, 90.317, 90.560)
      Space group: P 1
      U matrix:  {{-0.2595,  0.3410,  0.9035},
                  { 0.3839,  0.8949, -0.2275},
                  {-0.8862,  0.2878, -0.3632}}
      B matrix:  {{ 0.0172,  0.0000,  0.0000},
                  { 0.0002,  0.0171,  0.0000},
                  { 0.0001,  0.0001,  0.0067}}
      A = UB:    {{-0.0043,  0.0059,  0.0060},
                  { 0.0067,  0.0153, -0.0015},
                  {-0.0152,  0.0049, -0.0024}}


  757 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 1)
  ################################################################################


  Summary statistics for observations matched to predictions:
  -----------------------------------------------------------------------
  |                   | Min     | Q1       | Med      | Q3     | Max    |
  -----------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.7665 | -0.4922  | -0.05848 | 0.1489 | 0.4568 |
  | Yc - Yo (mm)      | -0.8621 | -0.4161  | 0.04831  | 0.2403 | 0.5781 |
  | Phic - Phio (deg) | -0.442  | -0.01297 | 0.1146   | 0.2693 | 0.9865 |
  | X weights         | 113.8   | 134.7    | 135      | 135.1  | 135.2  |
  | Y weights         | 119.2   | 134.9    | 135.1    | 135.2  | 135.2  |
  | Phi weights       | 162.5   | 177.1    | 177.5    | 177.7  | 177.8  |
  -----------------------------------------------------------------------


  Running refinement
  ------------------
  0 1 2 3 4

  Refinement steps
  ----------------
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.38369  | 0.37431  | 0.23548  |
  | 1    | 4049 | 0.12009  | 0.11387  | 0.18697  |
  | 2    | 4049 | 0.088057 | 0.081596 | 0.14271  |
  | 3    | 4049 | 0.048008 | 0.048841 | 0.076388 |
  | 4    | 4049 | 0.026475 | 0.035665 | 0.02821  |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.13752 | 0.18741 | 0.11533  |
  ---------------------------------------------
  Increasing resolution to 3.5 Angstrom
  model 1 (18471 reflections):
  Crystal:
      Unit cell: (57.783, 57.795, 149.981, 90.039, 90.024, 90.007)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3448,  0.9021},
                  { 0.3909,  0.8916, -0.2285},
                  {-0.8832,  0.2934, -0.3660}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  { 0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}


  86 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 2)
  ################################################################################


  Summary statistics for observations matched to predictions:
  -------------------------------------------------------------------------
  |                   | Min     | Q1       | Med       | Q3      | Max    |
  -------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.2894 | -0.04219 | -0.007729 | 0.01458 | 0.2122 |
  | Yc - Yo (mm)      | -0.7472 | -0.03726 | -0.0137   | 0.01058 | 0.2652 |
  | Phic - Phio (deg) | -1.045  | -0.01004 | 0.001227  | 0.01285 | 0.9063 |
  | X weights         | 110.6   | 134.7    | 135       | 135.1   | 135.2  |
  | Y weights         | 114     | 134.8    | 135.1     | 135.2   | 135.2  |
  | Phi weights       | 160.2   | 177.2    | 177.5     | 177.7   | 177.8  |
  -------------------------------------------------------------------------


  Running refinement
  ------------------
  0

  Refinement steps
  ----------------
  -----------------------------------------------
  | Step | Nref | RMSD_X  | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)    | (mm)     | (deg)    |
  -----------------------------------------------
  | 0    | 4049 | 0.04984 | 0.043948 | 0.02381  |
  -----------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.23649 | 0.23959 | 0.16135  |
  ---------------------------------------------
  Increasing resolution to 2.5 Angstrom
  model 1 (47547 reflections):
  Crystal:
      Unit cell: (57.767, 57.785, 149.992, 90.051, 90.004, 89.998)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3449,  0.9021},
                  { 0.3909,  0.8916, -0.2285},
                  {-0.8832,  0.2934, -0.3660}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}


  137 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 3)
  ################################################################################


  Summary statistics for observations matched to predictions:
  ------------------------------------------------------------------------
  |                   | Min     | Q1       | Med      | Q3      | Max    |
  ------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.4052 | -0.02696 | 0.008701 | 0.05057 | 0.3013 |
  | Yc - Yo (mm)      | -0.73   | -0.06316 | -0.02221 | 0.01147 | 0.2737 |
  | Phic - Phio (deg) | -1.039  | -0.01125 | 0.001069 | 0.01416 | 0.911  |
  | X weights         | 101.4   | 134.1    | 134.8    | 135.1   | 135.2  |
  | Y weights         | 103.4   | 134      | 134.8    | 135.1   | 135.2  |
  | Phi weights       | 157.8   | 176.8    | 177.4    | 177.7   | 177.8  |
  ------------------------------------------------------------------------


  Running refinement
  ------------------
  0 1 2 3 4

  Refinement steps
  ----------------
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.061077 | 0.064421 | 0.022755 |
  | 1    | 4049 | 0.05994  | 0.05691  | 0.022316 |
  | 2    | 4049 | 0.059291 | 0.055404 | 0.021625 |
  | 3    | 4049 | 0.057855 | 0.053071 | 0.020781 |
  | 4    | 4049 | 0.054599 | 0.049016 | 0.020268 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.28272 | 0.24168 | 0.13326  |
  ---------------------------------------------
  Increasing resolution to 1.5 Angstrom
  model 1 (113986 reflections):
  Crystal:
      Unit cell: (57.793, 57.802, 150.030, 90.024, 90.011, 89.995)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3451,  0.9020},
                  { 0.3910,  0.8915, -0.2287},
                  {-0.8831,  0.2934, -0.3661}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}


  328 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 4)
  ################################################################################


  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.4434 | -0.03989 | 0.001557   | 0.0502  | 0.5943 |
  | Yc - Yo (mm)      | -1.172  | -0.07835 | -0.02655   | 0.0125  | 1.477  |
  | Phic - Phio (deg) | -1.424  | -0.01516 | -0.0005431 | 0.01418 | 0.9079 |
  | X weights         | 81.12   | 131.3    | 133.8      | 134.9   | 135.2  |
  | Y weights         | 87.23   | 130      | 133.3      | 134.7   | 135.2  |
  | Phi weights       | 145.2   | 176.2    | 177.4      | 177.8   | 177.8  |
  --------------------------------------------------------------------------


  Running refinement
  ------------------
  0 1 2 3 4 5

  Refinement steps
  ----------------
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.077047 | 0.088022 | 0.02738  |
  | 1    | 4049 | 0.07408  | 0.075985 | 0.027343 |
  | 2    | 4049 | 0.072926 | 0.074442 | 0.027333 |
  | 3    | 4049 | 0.070029 | 0.070914 | 0.027238 |
  | 4    | 4049 | 0.06367  | 0.063312 | 0.027096 |
  | 5    | 4049 | 0.054825 | 0.052748 | 0.026893 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  --------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y | RMSD_Z   |
  |     |      | (px)    | (px)   | (images) |
  --------------------------------------------
  | 0   | 4049 | 0.28944 | 0.2715 | 0.17801  |
  --------------------------------------------
  Increasing resolution to 0.5 Angstrom
  model 1 (114691 reflections):
  Crystal:
      Unit cell: (57.785, 57.797, 150.018, 90.019, 90.000, 89.993)
      Space group: P 1
      U matrix:  {{-0.2591,  0.3453,  0.9020},
                  { 0.3911,  0.8914, -0.2290},
                  {-0.8831,  0.2934, -0.3661}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}


  341 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 5)
  ################################################################################


  Summary statistics for observations matched to predictions:
  ------------------------------------------------------------------------
  |                   | Min    | Q1       | Med       | Q3      | Max    |
  ------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.574 | -0.03382 | -0.004531 | 0.03071 | 0.6556 |
  | Yc - Yo (mm)      | -1.409 | -0.02517 | 0.003335  | 0.02913 | 1.261  |
  | Phic - Phio (deg) | -1.427 | -0.01379 | 0.0002674 | 0.01498 | 0.9102 |
  | X weights         | 81.12  | 131.2    | 133.8     | 134.9   | 135.2  |
  | Y weights         | 87.23  | 130      | 133.3     | 134.7   | 135.2  |
  | Phi weights       | 145.2  | 176.2    | 177.5     | 177.8   | 177.8  |
  ------------------------------------------------------------------------


  Running refinement
  ------------------
  0

  Refinement steps
  ----------------
  ------------------------------------------------
  | Step | Nref | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |      | (mm)     | (mm)     | (deg)    |
  ------------------------------------------------
  | 0    | 4049 | 0.050185 | 0.045811 | 0.028146 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 4049 | 0.29133 | 0.26534 | 0.18682  |
  ---------------------------------------------
  Final refined crystal models:
  model 1 (114691 reflections):
  Crystal:
      Unit cell: (57.782, 57.797, 150.019, 90.017, 90.000, 89.996)
      Space group: P 1
      U matrix:  {{-0.2591,  0.3454,  0.9020},
                  { 0.3911,  0.8914, -0.2290},
                  {-0.8831,  0.2934, -0.3661}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}



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
         9  0.0195 dg 0.070   4049   tetragonal tP  57.79  57.79 150.01  90.00  90.00  90.00  501004      a,b,c
         8  0.0195 dg 0.069   4049 orthorhombic oC  81.72  81.74 150.02  90.00  90.00  90.00 1002112  a-b,a+b,c
         7  0.0175 dg 0.069   4049 orthorhombic oP  57.78  57.80 150.01  90.00  90.00  90.00  500989      a,b,c
         6  0.0195 dg 0.068   4049   monoclinic mC  81.72  81.73 150.02  90.00  89.99  90.00 1002056  a-b,a+b,c
         5  0.0191 dg 0.069   4049   monoclinic mC  81.74  81.72 150.02  90.00  90.01  90.00 1002108 a+b,-a+b,c
         4  0.0175 dg 0.069   4049   monoclinic mP  57.78  57.80 150.01  90.00  90.00  90.00  500989      a,b,c
         3  0.0170 dg 0.069   4049   monoclinic mP  57.78 150.01  57.80  90.00  89.99  90.00  501034   -a,-c,-b
         2  0.0044 dg 0.068   4049   monoclinic mP  57.80  57.78 150.02  90.00  90.02  90.00  500980   -b,-a,-c
         1  0.0000 dg 0.068   4049    triclinic aP  57.78  57.80 150.02  90.02  90.00  90.00  501002      a,b,c
  -------------------------------------------------------------------------------------------------------------
  usr+sys time: 2.06 seconds, ticks: 3148462, micro-seconds/tick: 0.654

In this example we would continue processing (i.e. proceed to the refinement
step, perhaps) with :samp:`bravais_setting_9.json`. Sometimes it may be
necessary to reindex the :samp:`indexed.pickle` file output by dials.index.
However, in this case as the change of basis operator to the chosen setting
is the identity operator (:samp:`a,b,c`) this step is not needed::

  dials.reindex indexed.pickle change_of_basis_op=a,b,c

This outputs the file :samp:`reflections_reindexed.pickle` which should be
used as input to downstream programs in place of :samp:`indexed.pickle`.


Refinement
^^^^^^^^^^

Although the model is already refined in indexing we can also add a refinement
step in here to allow e.g. scan varying refinement as here.

::

  dials.refine bravais_setting_9.json reflections_reindexed.pickle \
  refinement.parameterisation.crystal.scan_varying=true \
  refinement.reflections.use_all_reflections=true

This one on the other hand would probably stand to be *more* verbose!

::

  Configuring refiner
  Performing refinement
  Saving refined experiments to refined_experiments.json

Integration
^^^^^^^^^^^

After the refinement is done the next step is integration, which is performed
by the program :doc:`dials.integrate </programs/dials_integrate>`.

::

  dials.integrate refined_experiments.json indexed.pickle

This program outputs a lot of information as integration progresses,
concluding with a summary of the integration results.

::

  Summary of integration results binned by resolution
   ---------------------------------------------------------------------------------------
   d min |  d max | # full | # part | # over | # ice | # sum | # prf | <I/sigI> | <I/sigI>
         |        |        |        |        |       |       |       |    (sum) |    (prf)
   ---------------------------------------------------------------------------------------
    1.17 |   1.19 |    304 |      2 |      0 |     0 |   306 |   233 |      0.4 |      0.5
    1.19 |   1.21 |   1063 |      5 |      0 |     0 |  1068 |   911 |      0.4 |      0.5
    1.21 |   1.23 |   2265 |     13 |      0 |     0 |  2278 |  2061 |      0.5 |      0.6
    1.23 |   1.26 |   3719 |     21 |      0 |     0 |  3740 |  3519 |      0.5 |      0.7
    1.26 |   1.28 |   5346 |     29 |      0 |     0 |  5375 |  5080 |      0.6 |      0.8
    1.28 |   1.31 |   7116 |     45 |      0 |     0 |  7161 |  6818 |      0.7 |      0.8
    1.31 |   1.35 |   9374 |     58 |      0 |     0 |  9432 |  9044 |      0.8 |      1.0
    1.35 |   1.38 |  12321 |     78 |      0 |     0 | 12399 | 11981 |      0.9 |      1.1
    1.38 |   1.42 |  16781 |     94 |      0 |     0 | 16875 | 16329 |      1.0 |      1.2
    1.42 |   1.47 |  19951 |    133 |      0 |     0 | 20084 | 19691 |      1.2 |      1.5
    1.47 |   1.52 |  23395 |    193 |      0 |     0 | 23588 | 23218 |      1.5 |      1.8
    1.52 |   1.58 |  23905 |    204 |      0 |     0 | 24109 | 23972 |      1.8 |      2.1
    1.58 |   1.66 |  25334 |    210 |      0 |     0 | 25544 | 25436 |      2.2 |      2.6
    1.66 |   1.74 |  24067 |    185 |      0 |     0 | 24252 | 24185 |      2.7 |      3.1
    1.74 |   1.85 |  24607 |    189 |      0 |     0 | 24796 | 24750 |      3.5 |      4.0
    1.85 |   2.00 |  25544 |    222 |      0 |     0 | 25766 | 25742 |      4.9 |      5.4
    2.00 |   2.20 |  24539 |    205 |      0 |     0 | 24744 | 24734 |      6.6 |      7.2
    2.20 |   2.52 |  25538 |    195 |      0 |     0 | 25733 | 25728 |      8.8 |      9.4
    2.52 |   3.17 |  25050 |    236 |      0 |     0 | 25286 | 25283 |     12.7 |     13.3
    3.17 | 151.26 |  25609 |    225 |      0 |     0 | 25834 | 25831 |     25.4 |     25.6
   ---------------------------------------------------------------------------------------

   Summary of integration results for the whole dataset
   ----------------------------------------------
   Number fully recorded                 | 370431
   Number partially recorded             | 4294
   Number with overloaded pixels         | 0
   Number in powder rings                | 0
   Number processed with summation       | 328370
   Number processed with profile fitting | 324546
   <I/sigI> (summation)                  | 5.6
   <I/sigI> (profile fitting)            | 6.1
   ----------------------------------------------

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

  Removing 22782 reflections with negative variance
  Removing 27509 profile reflections with negative variance
  Removing 1963 incomplete reflections
  Title: from dials.export_mtz
  Space group symbol from file: P4
  Space group number from file: 75
  Space group from matrices: P 4 (No. 75)
  Point group symbol from file: 4
  Number of batches: 540
  Number of crystals: 1
  Number of Miller indices: 322471
  Resolution range: 150.008 1.17004
  History:
  Crystal 1:
    Name: XTAL
    Project: DIALS
    Id: 1
    Unit cell: (57.7876, 57.7876, 150.008, 90, 90, 90)
    Number of datasets: 1
    Dataset 1:
      Name: FROMDIALS
      Id: 1
      Wavelength: 0.97625
      Number of columns: 14
      label        #valid  %valid    min     max type
      H            322471 100.00%   0.00   47.00 H: index h,k,l
      K            322471 100.00%   0.00   46.00 H: index h,k,l
      L            322471 100.00%   0.00  114.00 H: index h,k,l
      M_ISYM       322471 100.00%   1.00    8.00 Y: M/ISYM, packed partial/reject flag and symmetry number
      BATCH        322471 100.00%   2.00  539.00 B: BATCH number
      IPR          322471 100.00%  -2.51 2938.60 J: intensity
      SIGIPR       322471 100.00%   0.05   54.24 Q: standard deviation
      I            322471 100.00% -24.60 3059.63 J: intensity
      SIGI         322471 100.00%   0.09   55.45 Q: standard deviation
      FRACTIONCALC 322471 100.00%   1.00    1.00 R: real
      XDET         322471 100.00%   6.41 2456.14 R: real
      YDET         322471 100.00%   5.70 2520.45 R: real
      ROT          322471 100.00%  82.00  162.70 R: real
      LP           322471 100.00%   0.00    0.76 R: real

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

to get merged data for downstream analysis. The output from this will include
the merging statistics which will give some idea of the data quality. Often
passing in a sensible resolution limit to aimless is also helpful... this should
give you something like::

  Summary data for        Project: DIALS Crystal: XTAL Dataset: FROMDIALS

                                             Overall  InnerShell  OuterShell
  Low resolution limit                      150.01    150.01      1.32
  High resolution limit                       1.30      7.12      1.30

  Rmerge  (within I+/I-)                     0.064     0.024     0.405
  Rmerge  (all I+ and I-)                    0.072     0.026     0.474
  Rmeas (within I+/I-)                       0.079     0.030     0.559
  Rmeas (all I+ & I-)                        0.080     0.030     0.593
  Rpim (within I+/I-)                        0.045     0.017     0.384
  Rpim (all I+ & I-)                         0.034     0.014     0.349
  Rmerge in top intensity bin                0.029        -         -
  Total number of observations              307413      2248      5428
  Total number unique                        62326       499      2461
  Mean((I)/sd(I))                             10.8      26.6       1.5
  Mn(I) half-set correlation CC(1/2)         0.999     0.999     0.707
  Completeness                                98.2      99.8      79.6
  Multiplicity                                 4.9       4.5       2.2

  Anomalous completeness                      92.2     100.0      47.2
  Anomalous multiplicity                       2.4       3.0       1.5
  DelAnom correlation between half-sets      0.001     0.346    -0.042
  Mid-Slope of Anom Normal Probability       0.944       -         -

.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html
