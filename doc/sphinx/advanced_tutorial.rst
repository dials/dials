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
``/here/are/all/images*.cbf`` will do sensible processing, with a static model
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
      nproc = 4
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

  refinement {
    reflections {
      use_all_reflections = true
    }
  }
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
  | 0    | 7308 | 0.37964  | 0.37554  | 0.24355  |
  | 1    | 7308 | 0.12459  | 0.11211  | 0.19516  |
  | 2    | 7308 | 0.090179 | 0.079205 | 0.14681  |
  | 3    | 7308 | 0.048037 | 0.047679 | 0.076444 |
  | 4    | 7308 | 0.026579 | 0.036252 | 0.028052 |
  ------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |      | (px)    | (px)    | (images) |
  ---------------------------------------------
  | 0   | 7308 | 0.13967 | 0.19176 | 0.11965  |
  ---------------------------------------------
  Increasing resolution to 3.5 Angstrom
  model 1 (18473 reflections):
  Crystal:
      Unit cell: (57.784, 57.797, 149.985, 90.036, 90.022, 90.007)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3449,  0.9021},
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
  --------------------------------------------------------------------------
  |                   | Min     | Q1       | Med       | Q3       | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.29   | -0.04199 | -0.007764 | 0.01469  | 0.2129 |
  | Yc - Yo (mm)      | -0.7444 | -0.03736 | -0.01427  | 0.009862 | 0.2624 |
  | Phic - Phio (deg) | -1.047  | -0.01004 | 0.001048  | 0.01251  | 0.906  |
  | X weights         | 110.6   | 134.7    | 135       | 135.1    | 135.2  |
  | Y weights         | 114     | 134.8    | 135.1     | 135.2    | 135.2  |
  | Phi weights       | 160.2   | 177.2    | 177.5     | 177.7    | 177.8  |
  --------------------------------------------------------------------------


  Running refinement
  ------------------
  0

  Refinement steps
  ----------------
  -------------------------------------------------
  | Step | Nref  | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |       | (mm)     | (mm)     | (deg)    |
  -------------------------------------------------
  | 0    | 17725 | 0.049448 | 0.045068 | 0.021263 |
  -------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ----------------------------------------------
  | Exp | Nref  | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |       | (px)    | (px)    | (images) |
  ----------------------------------------------
  | 0   | 17725 | 0.23652 | 0.24496 | 0.14492  |
  ----------------------------------------------
  Increasing resolution to 2.5 Angstrom
  model 1 (47548 reflections):
  Crystal:
      Unit cell: (57.768, 57.784, 149.993, 90.048, 90.004, 90.000)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3449,  0.9021},
                  { 0.3909,  0.8916, -0.2285},
                  {-0.8832,  0.2934, -0.3660}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  { 0.0000,  0.0173,  0.0000},
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
  | Xc - Xo (mm)      | -0.4088 | -0.02791 | 0.007601 | 0.04941 | 0.2994 |
  | Yc - Yo (mm)      | -0.7251 | -0.06249 | -0.0221  | 0.01166 | 0.2708 |
  | Phic - Phio (deg) | -1.041  | -0.01083 | 0.001335 | 0.01398 | 0.9119 |
  | X weights         | 101.4   | 134.1    | 134.8    | 135.1   | 135.2  |
  | Y weights         | 103.4   | 134      | 134.8    | 135.1   | 135.2  |
  | Phi weights       | 157.8   | 176.8    | 177.4    | 177.7   | 177.8  |
  ------------------------------------------------------------------------


  Running refinement
  ------------------
  0 1 2 3 4

  Refinement steps
  ----------------
  -------------------------------------------------
  | Step | Nref  | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |       | (mm)     | (mm)     | (deg)    |
  -------------------------------------------------
  | 0    | 46528 | 0.061731 | 0.063425 | 0.02268  |
  | 1    | 46528 | 0.060736 | 0.056615 | 0.02223  |
  | 2    | 46528 | 0.060134 | 0.055245 | 0.021665 |
  | 3    | 46528 | 0.05872  | 0.053004 | 0.021021 |
  | 4    | 46528 | 0.055422 | 0.048926 | 0.020678 |
  -------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  ---------------------------------------------
  | Exp | Nref  | RMSD_X | RMSD_Y  | RMSD_Z   |
  |     |       | (px)   | (px)    | (images) |
  ---------------------------------------------
  | 0   | 46528 | 0.287  | 0.24118 | 0.1366   |
  ---------------------------------------------
  Increasing resolution to 1.5 Angstrom
  model 1 (113985 reflections):
  Crystal:
      Unit cell: (57.790, 57.804, 150.026, 90.022, 90.007, 89.995)
      Space group: P 1
      U matrix:  {{-0.2593,  0.3452,  0.9020},
                  { 0.3910,  0.8915, -0.2287},
                  {-0.8831,  0.2934, -0.3661}}
      B matrix:  {{ 0.0173,  0.0000,  0.0000},
                  {-0.0000,  0.0173,  0.0000},
                  { 0.0000,  0.0000,  0.0067}}
      A = UB:    {{-0.0045,  0.0060,  0.0060},
                  { 0.0068,  0.0154, -0.0015},
                  {-0.0153,  0.0051, -0.0024}}


  329 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 4)
  ################################################################################


  Summary statistics for observations matched to predictions:
  -------------------------------------------------------------------------
  |                   | Min    | Q1       | Med        | Q3      | Max    |
  -------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.443 | -0.04157 | -0.0004459 | 0.04768 | 0.5899 |
  | Yc - Yo (mm)      | -1.169 | -0.07645 | -0.02574   | 0.01291 | 1.269  |
  | Phic - Phio (deg) | -1.43  | -0.01476 | -0.0003432 | 0.01429 | 0.9075 |
  | X weights         | 81.12  | 131.3    | 133.8      | 134.9   | 135.2  |
  | Y weights         | 87.23  | 130      | 133.3      | 134.7   | 135.2  |
  | Phi weights       | 145.2  | 176.2    | 177.4      | 177.8   | 177.8  |
  -------------------------------------------------------------------------


  Running refinement
  ------------------
  0 1 2 3 4 5

  Refinement steps
  ----------------
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 112548 | 0.07768  | 0.088135 | 0.028143 |
  | 1    | 112548 | 0.074814 | 0.076596 | 0.028317 |
  | 2    | 112548 | 0.073593 | 0.0749   | 0.02827  |
  | 3    | 112548 | 0.070528 | 0.071118 | 0.028156 |
  | 4    | 112548 | 0.063893 | 0.063233 | 0.028053 |
  | 5    | 112548 | 0.054902 | 0.052613 | 0.027942 |
  --------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 112548 | 0.29089 | 0.27229 | 0.18559  |
  -----------------------------------------------
  Increasing resolution to 0.5 Angstrom
  model 1 (114691 reflections):
  Crystal:
      Unit cell: (57.779, 57.795, 150.013, 90.016, 90.002, 89.995)
      Space group: P 1
      U matrix:  {{-0.2591,  0.3454,  0.9020},
                  { 0.3911,  0.8914, -0.2290},
                  {-0.8831,  0.2934, -0.3660}}
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
  -------------------------------------------------------------------------
  |                   | Min     | Q1       | Med        | Q3      | Max   |
  -------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.5625 | -0.03212 | -0.002373  | 0.03215 | 0.653 |
  | Yc - Yo (mm)      | -1.408  | -0.02572 | 0.002198   | 0.02784 | 1.249 |
  | Phic - Phio (deg) | -1.417  | -0.01396 | -2.606e-05 | 0.01474 | 0.908 |
  | X weights         | 81.12   | 131.2    | 133.8      | 134.9   | 135.2 |
  | Y weights         | 87.23   | 130      | 133.3      | 134.7   | 135.2 |
  | Phi weights       | 145.2   | 176.2    | 177.5      | 177.8   | 177.8 |
  -------------------------------------------------------------------------


  Running refinement
  ------------------
  0

  Refinement steps
  ----------------
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 113249 | 0.050333 | 0.047544 | 0.027976 |
  --------------------------------------------------
  RMSD target achieved

  RMSDs by experiment
  -------------------
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 113249 | 0.29253 | 0.27634 | 0.1865   |
  -----------------------------------------------
  Final refined crystal models:
  model 1 (114691 reflections):
  Crystal:
      Unit cell: (57.779, 57.795, 150.012, 90.016, 90.001, 89.995)
      Space group: P 1
      U matrix:  {{-0.2591,  0.3454,  0.9020},
                  { 0.3911,  0.8914, -0.2290},
                  {-0.8831,  0.2934, -0.3660}}
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
         9  0.0197 dg 0.069   4049   tetragonal tP  57.79  57.79 150.01  90.00  90.00  90.00  500936      a,b,c
         8  0.0197 dg 0.069   4049 orthorhombic oC  81.72  81.73 150.01  90.00  90.00  90.00 1001961  a-b,a+b,c
         7  0.0167 dg 0.069   4049 orthorhombic oP  57.78  57.79 150.01  90.00  90.00  90.00  500920      a,b,c
         6  0.0197 dg 0.068   4049   monoclinic mC  81.72  81.73 150.01  90.00  89.99  90.00 1001915  a-b,a+b,c
         5  0.0184 dg 0.069   4049   monoclinic mC  81.73  81.72 150.01  90.00  90.01  90.00 1001960 a+b,-a+b,c
         4  0.0167 dg 0.069   4049   monoclinic mP  57.78  57.79 150.01  90.00  90.00  90.00  500920      a,b,c
         3  0.0160 dg 0.069   4049   monoclinic mP  57.78 150.01  57.80  90.00  89.99  90.00  500960   -a,-c,-b
         2  0.0051 dg 0.067   4049   monoclinic mP  57.79  57.78 150.01  90.00  90.01  90.00  500911   -b,-a,-c
         1  0.0000 dg 0.067   4049    triclinic aP  57.78  57.79 150.01  90.01  90.00  90.00  500930      a,b,c
  -------------------------------------------------------------------------------------------------------------
  usr+sys time: 0.84 seconds
  wall clock time: 3.92 seconds

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
step using :doc:`dials.refine </programs/dials_refine>` in here to allow e.g.
scan-varying refinement. In fact, this
dataset is of such good quality that with default options scan-varying
refinement would terminate immediately because the RMSD target was already
achieved during the indexing step. To ensure we squeeze the best possible results
from this data we use the expert parameter ``bin_size_fraction`` to set the
RMSD target to zero in each dimension. This ensures that refinement continues
until RMSD convergence.

As an aside, to show all the options up to and including ``expert_level = 1``
use this command::

  dials.refine -c -e 1

Now, our refinement job is specified as::

  dials.refine bravais_setting_9.json reflections_reindexed.pickle \
  refinement.parameterisation.crystal.scan_varying=true \
  refinement.reflections.use_all_reflections=true \
  refinement.target.bin_size_fraction=0.0

The main product of this is the file ``refined_experiments.json``

::

  The following parameters have been modified:

  refinement {
    parameterisation {
      crystal {
        scan_varying = true
      }
    }
    target {
      bin_size_fraction = 0.0
    }
    reflections {
      use_all_reflections = true
    }
  }
  input {
    experiments = bravais_setting_9.json
    reflections = reflections_reindexed.pickle
  }

  Configuring refiner

  Summary statistics for observations matched to predictions:
  --------------------------------------------------------------------------
  |                   | Min    | Q1        | Med        | Q3      | Max    |
  --------------------------------------------------------------------------
  | Xc - Xo (mm)      | -0.55  | -0.03375  | -0.00355   | 0.0305  | 0.6397 |
  | Yc - Yo (mm)      | -1.405 | -0.03049  | -0.0006985 | 0.02968 | 1.227  |
  | Phic - Phio (deg) | -1.357 | -0.007731 | 0.007483   | 0.02335 | 0.9143 |
  | X weights         | 81.12  | 131.2     | 133.8      | 134.9   | 135.2  |
  | Y weights         | 87.23  | 130       | 133.3      | 134.7   | 135.2  |
  | Phi weights       | 145.2  | 176.2     | 177.5      | 177.8   | 177.8  |
  --------------------------------------------------------------------------

  Performing refinement

  Running refinement
  ------------------
  0 1 2 3 4 5 6 7 8

  Refinement steps
  ----------------
  --------------------------------------------------
  | Step | Nref   | RMSD_X   | RMSD_Y   | RMSD_Phi |
  |      |        | (mm)     | (mm)     | (deg)    |
  --------------------------------------------------
  | 0    | 113249 | 0.050561 | 0.050785 | 0.030271 |
  | 1    | 113249 | 0.05059  | 0.047453 | 0.029408 |
  | 2    | 113249 | 0.050682 | 0.046859 | 0.028885 |
  | 3    | 113249 | 0.050534 | 0.046412 | 0.0284   |
  | 4    | 113249 | 0.050215 | 0.046259 | 0.028003 |
  | 5    | 113249 | 0.050018 | 0.046328 | 0.027803 |
  | 6    | 113249 | 0.049967 | 0.046377 | 0.027738 |
  | 7    | 113249 | 0.049961 | 0.046383 | 0.027729 |
  | 8    | 113249 | 0.049961 | 0.046383 | 0.027728 |
  --------------------------------------------------
  RMSD no longer decreasing

  RMSDs by experiment
  -------------------
  -----------------------------------------------
  | Exp | Nref   | RMSD_X  | RMSD_Y  | RMSD_Z   |
  |     |        | (px)    | (px)    | (images) |
  -----------------------------------------------
  | 0   | 113249 | 0.29046 | 0.26967 | 0.18484  |
  -----------------------------------------------
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
  ----------------------------------------------------------------------------------------------------------
  d min |  d max | # full | # part | # over | # ice | # sum | # prf | <Ibg> | <I/sigI> | <I/sigI> | <CC prf>
        |        |        |        |        |       |       |       |       |    (sum) |    (prf) |
  ----------------------------------------------------------------------------------------------------------
  1.17 |   1.19 |    296 |      2 |      0 |     0 |   298 |   232 |  0.00 |     3.24 |     2.15 |     0.13
  1.19 |   1.21 |   1062 |      5 |      0 |     0 |  1067 |   935 |  0.00 |     3.30 |     2.20 |     0.11
  1.21 |   1.23 |   2265 |     13 |      0 |     0 |  2278 |  2086 |  0.00 |     3.38 |     2.25 |     0.12
  1.23 |   1.26 |   3713 |     20 |      0 |     0 |  3733 |  3553 |  0.00 |     3.43 |     2.32 |     0.14
  1.26 |   1.28 |   5340 |     30 |      0 |     0 |  5370 |  5143 |  0.00 |     3.49 |     2.39 |     0.16
  1.28 |   1.31 |   7113 |     44 |      0 |     0 |  7157 |  6914 |  0.00 |     3.54 |     2.43 |     0.17
  1.31 |   1.35 |   9365 |     56 |      0 |     0 |  9421 |  9168 |  0.00 |     3.66 |     2.54 |     0.20
  1.35 |   1.38 |  12326 |     77 |      0 |     0 | 12403 | 12097 |  0.00 |     3.69 |     2.63 |     0.24
  1.38 |   1.42 |  16762 |     95 |      0 |     0 | 16857 | 16496 |  0.01 |     3.53 |     2.58 |     0.25
  1.42 |   1.47 |  19946 |    137 |      0 |     0 | 20083 | 19856 |  0.02 |     3.48 |     2.65 |     0.29
  1.47 |   1.52 |  23305 |    459 |      0 |     0 | 23764 | 23511 |  0.03 |     3.41 |     2.73 |     0.33
  1.52 |   1.58 |  23792 |    539 |      0 |     0 | 24331 | 24288 |  0.05 |     3.30 |     2.82 |     0.38
  1.58 |   1.66 |  25230 |    522 |      0 |     0 | 25752 | 25706 |  0.07 |     3.23 |     3.01 |     0.44
  1.66 |   1.74 |  23971 |    472 |      0 |     0 | 24443 | 24406 |  0.09 |     3.34 |     3.36 |     0.50
  1.74 |   1.85 |  24507 |    461 |      0 |     0 | 24968 | 24943 |  0.12 |     3.88 |     4.05 |     0.57
  1.85 |   2.00 |  25444 |    514 |      0 |     0 | 25958 | 25932 |  0.16 |     5.29 |     5.45 |     0.64
  2.00 |   2.20 |  24468 |    440 |      0 |     0 | 24908 | 24891 |  0.20 |     7.09 |     7.23 |     0.70
  2.20 |   2.52 |  25445 |    445 |      0 |     0 | 25890 | 25866 |  0.23 |     9.28 |     9.43 |     0.74
  2.52 |   3.17 |  24970 |    486 |      0 |     0 | 25456 | 25420 |  0.31 |    12.97 |    13.18 |     0.76
  3.17 | 151.26 |  25503 |    595 |      0 |     0 | 26098 | 26071 |  0.38 |    25.41 |    25.26 |     0.76
  ----------------------------------------------------------------------------------------------------------

  Summary of integration results for the whole dataset
  ----------------------------------------------
  Number fully recorded                 | 369259
  Number partially recorded             | 8630
  Number with overloaded pixels         | 0
  Number in powder rings                | 0
  Number processed with summation       | 330235
  Number processed with profile fitting | 327514
  <Ibg>                                 | 0.13
  <I/sigI> (summation)                  | 6.81
  <I/sigI> (profile fitting)            | 6.56
  <CC prf>                              | 0.44
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

  Removing 24192 reflections with negative variance
  Removing 26869 profile reflections with negative variance
  Removing 4152 incomplete reflections
  Title: from dials.export_mtz
  Space group symbol from file: P4
  Space group number from file: 75
  Space group from matrices: P 4 (No. 75)
  Point group symbol from file: 4
  Number of batches: 540
  Number of crystals: 1
  Number of Miller indices: 322676
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
      label        #valid  %valid   min     max type
      H            322676 100.00%  0.00   47.00 H: index h,k,l
      K            322676 100.00%  0.00   46.00 H: index h,k,l
      L            322676 100.00%  0.00  114.00 H: index h,k,l
      M_ISYM       322676 100.00%  1.00    8.00 Y: M/ISYM, packed partial/reject flag and symmetry number
      BATCH        322676 100.00%  2.00  539.00 B: BATCH number
      IPR          322676 100.00% -1.79 2892.40 J: intensity
      SIGIPR       322676 100.00%  0.00   53.81 Q: standard deviation
      I            322676 100.00% -7.09 3060.88 J: intensity
      SIGI         322676 100.00%  0.07   55.45 Q: standard deviation
      FRACTIONCALC 322676 100.00%  1.00    1.00 R: real
      XDET         322676 100.00%  6.50 2456.34 R: real
      YDET         322676 100.00%  5.79 2520.69 R: real
      ROT          322676 100.00% 82.01  162.70 R: real
      LP           322676 100.00%  0.00    0.76 R: real

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
  Low resolution limit                      149.99    149.99      1.32
  High resolution limit                       1.30      7.12      1.30

  Rmerge  (within I+/I-)                     0.065     0.024     0.223
  Rmerge  (all I+ and I-)                    0.074     0.026     0.278
  Rmeas (within I+/I-)                       0.080     0.029     0.307
  Rmeas (all I+ & I-)                        0.083     0.029     0.346
  Rpim (within I+/I-)                        0.046     0.016     0.211
  Rpim (all I+ & I-)                         0.035     0.013     0.202
  Rmerge in top intensity bin                0.029        -         -
  Total number of observations              307414      2242      5485
  Total number unique                        62340       499      2472
  Mean((I)/sd(I))                             11.8      27.3       4.3
  Mn(I) half-set correlation CC(1/2)         0.999     0.999     0.663
  Completeness                                98.2      99.8      80.1
  Multiplicity                                 4.9       4.5       2.2

  Anomalous completeness                      92.3     100.0      47.8
  Anomalous multiplicity                       2.4       3.0       1.5
  DelAnom correlation between half-sets     -0.007     0.193     0.033
  Mid-Slope of Anom Normal Probability       1.097       -         -

.. _pointless: http://www.ccp4.ac.uk/html/pointless.html
.. _aimless: http://www.ccp4.ac.uk/html/aimless.html
.. _ctruncate: http://www.ccp4.ac.uk/html/ctruncate.html
