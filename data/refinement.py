#!/usr/bin/env python
#
# refinement.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from libtbx.phil import parse

phil_scope = parse('''

refinement
  .help = "Parameters to configure the refinement"
{

  verbosity = 1
    .help = "verbosity level"
    .type = int(value_min=0)

  go_fast = True
    .help = "set False to revert to old classes for scan-static refinement"
            "that do a loop over all reflections in Python rather than the"
            "vectorised version using flex array operations"
    .type = bool

  parameterisation
    .help = "Parameters to control the parameterisation of experimental models"
  {
    beam
      .help = "beam parameters"
    {
      fix = all *in_spindle_plane out_spindle_plane
        .help = "Whether to fix beam parameters. By default,"
                "in_spindle_plane is selected, and one of the two"
                "parameters is fixed. If a goniometer is present"
                "this leads to the beam orientation being restricted"
                "to a direction in the initial spindle-beam plane"
        .type = choice

      fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"
    }

    crystal
      .help = "crystal parameters"
    {
      fix = all cell orientation
        .help = "Fix crystal parameters"
        .type = choice

      cell_fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"

      orientation_fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"

      scan_varying = False
        .help = "Parameterise the crystal to vary during the scan"
        .type = bool

      num_intervals = *fixed_width absolute
        .help = "Choose the way to determine the number of intervals for scan-"
                "varying refinement"
        .type = choice

      interval_width_degrees = 36.0
        .help = "Width of scan between checkpoints in degrees"
        .type = float(value_min=0.)

      absolute_num_intervals = 5
        .help = "Number of intervals between checkpoints if scan_varying"
                "refinement is requested"
        .type = int(value_min=1)
     }

    detector
      .help = "detector parameters"
    {
      panels = *automatic single multiple hierarchical
        .help = "Select appropriate detector parameterisation. Both the"
                "single and multiple panel detector options treat the whole"
                "detector as a rigid body. The hierarchical parameterisation"
                "treats groups of panels as separate rigid bodies."
        .type = choice
      hierarchy_level = 0
        .help = "Level of the detector hierarchy (starting from the root at 0)"
                "at which to determine panel groups to parameterise independently"
        .type = int(value_min=0)
      fix = all position orientation
        .help = "Fix detector parameters. The translational parameters"
                "(position) may be set separately to the orientation."
        .type = choice
      fix_list = None
        .type = ints(value_min=0)
        .help = "Fix specified parameters by a list of indices"
    }
  }

  refinery
    .help = "Parameters to configure the refinery"
  {
    engine = SimpleLBFGS LBFGScurvs GaussNewton *LevMar
      .help = "The minimisation engine to use"
      .type = choice

    track_step = False
      .help = "Record parameter shifts history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_gradient = False
      .help = "Record parameter gradients history in the refinement journal, if"
              "the engine supports it."
      .type = bool

    track_parameter_correlation = False
      .help = "Record correlation matrix between columns of the Jacobian for"
              "each step of refinement."
      .type = bool

    log = None
      .help = "Filename for an optional log that a minimisation engine may use"
              "to write additional information"
      .type = path

    max_iterations = None
      .help = "Maximum number of iterations in refinement before termination."
              "None implies the engine supplies its own default."
      .type = int(value_min=1)
  }

  target
    .help = "Parameters to configure the target function"
  {

    rmsd_cutoff = *fraction_of_bin_size absolute
      .help = "Method to choose rmsd cutoffs. This is currently either as a"
              "fraction of the discrete units of the spot positional data, i.e."
              "(pixel width, pixel height, image thickness in phi), or a tuple"
              "of absolute values to use as the cutoffs"
      .type = choice

    bin_size_fraction = 0.33333
      .help = "Cut off in the natural discrete units of positional data, viz.,"
              "(pixel width, pixel height, image thickness in phi) to use to"
              "determine when the RMSD target is achieved. Only used if"
              "rmsd_cutoff = fraction_of_bin_size."
      .type = float(value_min=0.)

    absolute_cutoffs = None
      .help = "Absolute Values for the RMSD target achieved cutoffs in X, Y and"
              "Phi. The units are (mm, mm, rad)."
      .type = floats(size=3, value_min=0.)
  }

  reflections
    .help = "Parameters used by the reflection manager"
  {

    reflections_per_degree = 50
      .help = "The number of centroids per degree of the sweep to use in"
              "refinement."
      .type = float(value_min=0.)

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input"
              "reflection list is done"
      .type = int

    maximum_sample_size = None
      .help = "The maximum number of reflections to use in refinement."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)

    use_all_reflections = False
      .help = "Override reflections_per_degree and use all available centroids"
              "in refinement."
      .type = bool

    random_seed = 42
      .help = "Random seed to use when sampling to create a working set of"
              "reflections. May be int or None."
      .type = int

    minimum_number_of_reflections = 20
      .help = "The minimum number of input observations to allow a reflection"
              "manager to be constructed for."
      .type = int(value_min=1)

    close_to_spindle_cutoff = 0.1
      .help = "The inclusion criterion currently uses the volume of the"
              "parallelepiped formed by the spindle axis, the incident"
              "beam and the scattered beam. If this is lower than some"
              "value then the reflection is excluded from refinement."
              "In detector space, these are the reflections located close"
              "to the rotation axis."
      .type = float(value_min = 0)

    do_outlier_rejection = False
      .help = "Whether refinement should attempt to reject outliers. Warning:"
              "by default we assume the user code will have already done this!"
      .type = bool

    iqr_multiplier = 1.5
      .help = "The IQR multiplier used to detect outliers. A value of 1.5 gives"
              "Tukey's rule for outlier detection"
      .type = float(value_min = 0.)
  }
}

''')
