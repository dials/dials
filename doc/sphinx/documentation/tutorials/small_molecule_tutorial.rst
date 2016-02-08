Running the Individual Steps: Small Molecule
--------------------------------------------

FIXME write this section too, and in each section detail the available
parameters. Guess it would also be rather nice to write about the available
parameters, and :samp:`dials.parameters`. For this example will process some
small-molecule data (sugar) to really demonstrate how the command-line options
work.

Import
^^^^^^

The first stage of step-by-step DIALS processing is to import the data - all
that happens here is that the image headers are read, and a file describing
their contents (:samp:`datablock.json`) is written. It's worth noting that if
this file is changed subsequent processing can use this.

::

  dials.import ~/data/sugar/15_n_120K_1_00*cbf

The output just describes what the software understands of the images it was
passed - not very interesting but useful to make sure it all makes sense.

::

  --------------------------------------------------------------------------------
  DataBlock 0
    format: <class 'dxtbx.format.FormatCBFFullPilatus.FormatCBFFullPilatus'>
    num images: 645
    num sweeps: 1
    num stills: 0
  --------------------------------------------------------------------------------
  Writing datablocks to datablock.json

Find Spots
^^^^^^^^^^

The first "real" task in any DIALS processing will be the spot finding - while
there are plenty of options the defaults often seem to do sensible things.

::

  dials.find_spots datablock.json

This will just report the number of spots found - guess we could probably write
some more interesting output.

::

  Configuring spot finder from input parameters
  --------------------------------------------------------------------------------
  Finding strong spots in imageset 0
  --------------------------------------------------------------------------------

  Finding spots in image 0 to 645...
  Extracted strong pixels from images.......................................24.35s
  Merged 8 pixel lists with 34901 pixels.....................................0.00s
  Extracted 6951 spots.......................................................0.09s
  Calculated 6951 spot centroids.............................................0.09s
  Calculated 6951 spot intensities...........................................0.01s
  Filtered 1448 spots by number of pixels....................................0.00s
  Filtered 1374 spots by peak-centroid distance..............................0.00s

  --------------------------------------------------------------------------------
  Saved 1374 reflections to strong.pickle....................................0.01s
  Total time:  26.5041749477

Indexing
^^^^^^^^

Here we will use a bunch more options to properly index the data, as the unit
cell from sucrose is very small! Method of fft1d corresponds to the 1D FFT
indexing rather than the 3D FFT default, the grid spacing clues just determine
how the indexing is done. N.B. the defaults are configured for macromolecular
crystallography.

::

  dials.index method=fft1d reciprocal_space_grid.d_min=1
    refinement_protocol.d_min_start=1 datablock.json strong.pickle max_cell=20

The output for this is rather verbose: FIXME perhaps I should abbreviate it
some?

::

  reference {
    detector = None
    beam = None
  }
  discover_better_experimental_model = False
  min_cell = 20
  max_cell = 20
  reciprocal_space_grid {
    n_points = 256
    d_min = 1
  }
  sigma_phi_deg = None
  b_iso = 200
  rmsd_cutoff = 15
  scan_range = None
  known_symmetry {
    space_group = None
    unit_cell = None
    relative_length_tolerance = 0.1
    absolute_angle_tolerance = 10
  }
  optimise_initial_basis_vectors = False
  debug = False
  debug_plots = False
  show_timing = False
  refinement {
    parameterisation {
      beam {
        fix = all *in_spindle_plane out_spindle_plane
        fix_list = None
      }
      crystal {
        fix = all cell orientation
        cell_fix_list = None
        orientation_fix_list = None
        scan_varying = False
        num_intervals = *fixed_width absolute
        interval_width_degrees = 36.0
        absolute_num_intervals = 5
      }
      detector {
        panels = *automatic single multiple hierarchical
        hierarchy_level = 0
        fix = all position orientation
        fix_list = None
      }
    }
    refinery {
      engine = SimpleLBFGS LBFGScurvs GaussNewtonIterations *LevMarIterations
      track_step = False
      track_gradient = False
      track_parameter_correlation = False
      log = None
      max_iterations = None
    }
    target {
      rmsd_cutoff = *fraction_of_bin_size absolute
      bin_size_fraction = 0.33333
      absolute_cutoffs = None
    }
    reflections {
      reflections_per_degree = 50
      minimum_sample_size = 1000
      maximum_number_of_reflections = None
      use_all_reflections = False
      random_seed = 42
      minimum_number_of_reflections = 20
      close_to_spindle_cutoff = 0.1
      do_outlier_rejection = False
      iqr_multiplier = 1.5
    }
  }
  refinement_protocol {
    weight_outlier_n_sigma = 5
    n_macro_cycles = 3
    d_min_step = 1.0
    d_min_start = 1
    d_min_final = None
    verbosity = 1
    outlier_rejection {
      hkl_tolerance = 0.3
    }
  }
  method = fft3d *fft1d real_space_grid_search
  multiple_lattice_search {
    cluster_analysis_search = False
    recycle_unindexed_reflections = False
    recycle_unindexed_reflections_cutoff = 0.1
    max_lattices = None
    cluster_analysis {
      method = *dbscan hcluster
      hcluster {
        linkage {
          method = *ward
          metric = *euclidean
        }
        cutoff = 15
        cutoff_criterion = *distance inconsistent
      }
      dbscan {
        eps = 0.05
        min_samples = 30
      }
      min_cluster_size = 20
      intersection_union_ratio_cutoff = 0.4
    }
  }
  Detector:
  Panel:
    pixel_size:{0.172,0.172}
    image_size: {487,619}
    trusted_range: {-1,243592}
    fast_axis: {1,0,0}
    slow_axis: {0,-0.866025,-0.5}
    origin: {-41.05,87.8104,-47.4522}

  Scan:
      image range:   {1,645}
      oscillation:   {-92,0.2}

  Goniometer:
      Rotation axis:  {1,-3.17778e-14,1.59583e-14}
      Fixed rotation: {0.661179,0.0297045,0.74964,-0.284305,-0.914768,0.287003,0.694272,-0.402887,-0.59638}

  Beam:
      wavelength: 0.6889
      sample to source direction : {0,0,1}
      divergence: 0
      sigma divergence: 0
      polarization normal: {0,1,0}
      polarization fraction: 0.8

  model 1 (136 reflections):
  Crystal:
      Unit cell: (7.493, 8.354, 10.337, 89.635, 77.719, 88.128)
      Space group: P 1
      U matrix:  {{ 0.5248, -0.3912,  0.7560},
                  {-0.2594, -0.9194, -0.2956},
                  { 0.8107, -0.0410, -0.5840}}
      B matrix:  {{ 0.1335,  0.0000,  0.0000},
                  {-0.0044,  0.1198,  0.0000},
                  {-0.0291,  0.0001,  0.0990}}
      A = UB:    {{ 0.0498, -0.0468,  0.0749},
                  {-0.0220, -0.1101, -0.0293},
                  { 0.1254, -0.0050, -0.0578}}


  801 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 1)
  ################################################################################


  Running refinement
  ------------------
  0 1 2 3 4 5 6 7 8 9 10 11

  Refinement steps
  ----------------
  Step Nref Objective RMSD_X RMSD_Y RMSD_Phi
  0 132 47422 0.4368 1.9379 0.024119
  1 132 6910.7 0.36193 0.3461 0.01491
  2 132 5063.1 0.34305 0.33002 0.012053
  3 132 3846.4 0.32599 0.32141 0.009753
  4 132 3396 0.31899 0.32687 0.0086288
  5 132 3313.8 0.31916 0.3268 0.0084046
  6 132 3293.2 0.32082 0.32416 0.0083628
  7 132 3270.3 0.32025 0.32331 0.008318
  8 132 3251.2 0.31897 0.32334 0.0082794
  9 132 3246.1 0.31816 0.32369 0.0082696
  10 132 3245.8 0.31794 0.32383 0.0082692
  11 132 3245.8 0.31792 0.32384 0.0082692
  RMSD no longer decreasing
  Increasing resolution to 0.0 Angstrom
  model 1 (1034 reflections):
  Crystal:
      Unit cell: (7.653, 8.586, 10.666, 89.631, 77.289, 89.852)
      Space group: P 1
      U matrix:  {{ 0.5443, -0.3908,  0.7423},
                  {-0.2534, -0.9201, -0.2986},
                  { 0.7997, -0.0255, -0.5999}}
      B matrix:  {{ 0.1307,  0.0000,  0.0000},
                  {-0.0003,  0.1165,  0.0000},
                  {-0.0295, -0.0007,  0.0961}}
      A = UB:    {{ 0.0494, -0.0460,  0.0713},
                  {-0.0240, -0.1070, -0.0287},
                  { 0.1222, -0.0026, -0.0577}}


  62 unindexed reflections

  ################################################################################
  Starting refinement (macro-cycle 2)
  ################################################################################


  Running refinement
  ------------------
  0 1 2 3 4 5 6 7 8 9 10

  Refinement steps
  ----------------
  Step Nref Objective RMSD_X RMSD_Y RMSD_Phi
  0 1027 50906 0.46797 0.43965 0.011741
  1 1027 31234 0.30073 0.31112 0.010643
  2 1027 27298 0.25668 0.27272 0.010466
  3 1027 23548 0.2175 0.21235 0.010339
  4 1027 21868 0.20436 0.17651 0.010241
  5 1027 21595 0.20583 0.17004 0.010195
  6 1027 21456 0.2056 0.16845 0.010167
  7 1027 21368 0.20484 0.16744 0.010154
  8 1027 21353 0.20449 0.16708 0.010154
  9 1027 21352 0.20442 0.16702 0.010155
  10 1027 21352 0.20442 0.16701 0.010155
  RMSD no longer decreasing
  Final refined crystal models:
  model 1 (1034 reflections):
  Crystal:
      Unit cell: (7.743, 8.707, 10.817, 90.180, 76.968, 90.281)
      Space group: P 1
      U matrix:  {{ 0.5416, -0.3876,  0.7460},
                  {-0.2449, -0.9216, -0.3011},
                  { 0.8042, -0.0196, -0.5940}}
      B matrix:  {{ 0.1291,  0.0000,  0.0000},
                  { 0.0006,  0.1148,  0.0000},
                  {-0.0299,  0.0002,  0.0949}}
      A = UB:    {{ 0.0474, -0.0443,  0.0708},
                  {-0.0232, -0.1059, -0.0286},
                  { 0.1216, -0.0024, -0.0564}}

  usr+sys time: 20.34 seconds, ticks: 99507393, micro-seconds/tick: 0.204
  wall clock time: 20.48 seconds

Refinement
^^^^^^^^^^

Although the model is already refined in indexing we can also add a refineent
step in here to allow e.g. scan varying refinement (though with this data we are
unlikely to really have enough measurements to do this!)

::

  dials.refine experiments.json indexed.pickle

This one on the other hand would probably stand to be *more* verbose!

::

  Configuring refiner
  Performing refinement
  Saving refined experiments to refined_experiments.json

Integration
^^^^^^^^^^^

After refinement is complete integration may proceed - though in this example it
failed so this is not a great example!

::

  %%%FIXME make command line right!

  dials.integrate refined_experiments.json

  FIXME ADD IN HERE RESULTS

Exporting as MTZ
^^^^^^^^^^^^^^^^

::

  FIXME ADD IN HERE EXPORT STUFF

