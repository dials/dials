"""
Phil scope of options for scaling.
"""
import iotbx.phil

phil_scope = iotbx.phil.parse('''
  overwrite_existing_models = False
    .type = bool
    .help = "If True, create new scaling models for all datasets"
    .expert_level = 0
  parameterisation {
    scale_term = True
      .type = bool
      .help = "Option to turn off scale correction (for physical/KB
               default models)."
      .expert_level = 2
    scale_interval = 15.0
      .type = float(value_min=1.0)
      .help = "Rotation (phi) interval between model parameters for the scale
               component (physical model)."
      .expert_level = 1
    decay_term = True
      .type = bool
      .help = "Option to turn off decay correction (for physical/array/KB
               default models)."
      .expert_level = 1
    decay_interval = 20.0
      .type = float(value_min=1.0)
      .help = "Rotation (phi) interval between model parameters for the decay
               component (physical/array default models)."
      .expert_level = 1
    n_resolution_bins = 10
      .type = int(value_min=1)
      .help = "Number of resolution bins to use for the decay term in the
               array-based model."
      .expert_level = 1
    decay_restraint = 1e-1
      .type = float(value_min=0.0)
      .help = "Weight to weakly restrain B-values to 0 for physical model."
      .expert_level = 2
    absorption_term = True
      .type = bool
      .help = "Option to turn off absorption correction (for physical/array
               default models)."
      .expert_level = 1
    lmax = 4
      .type = int(value_min=2)
      .help = "Number of spherical harmonics to include for absorption
              correction (for physical default model), recommended to be no
              more than 6."
      .expert_level = 1
    surface_weight = 1e6
      .type = float(value_min=0.0)
      .help = "Restraint weight applied to spherical harmonic terms in the
               physical model absorption correction."
      .expert_level = 1
    modulation_term = False
      .type = bool
      .help = "Option to turn on a detector correction for the array default
               model."
      .expert_level = 2
    n_modulation_bins = 20
      .type = int(value_min=1)
      .help = "Number of bins in each dimension (applied to both x and y) for
              binning the detector position for the modulation term of the
              array model."
      .expert_level = 2
    n_absorption_bins = 3
      .type = int(value_min=1)
      .help = "Number of bins in each dimension (applied to both x and y) for
              binning the detector position for the absorption term of the
              array model."
      .expert_level = 1
  }
  reflection_selection {
    E2_range = 0.8, 5.0
      .type = floats(size=2)
      .help = "Minimum and maximum normalised E^2 value to used to select a
              subset of reflections for minimising the scaling model."
      .expert_level = 1
    Isigma_range = -5.0, 0.0
      .type = floats(size=2)
      .help = "Minimum and maximum I/sigma values used to subset of reflections
              to determine the scaling model. Here a value of 0.0 for the max
              means no limit applied."
      .expert_level = 1
    d_range = None
      .type = floats(size=2)
      .help = "Minimum and maximum - values used to subset of reflections
              to determine the scaling model."
      .expert_level = 1
    min_partiality = 0.6
      .type = float
      .help = "Minimum partiality to use when selecting subset of reflections
               to determine the scaling model."
      .expert_level = 2
    intensity_choice = profile sum *combine
      .type = choice
      .help = "Option to choose from profile fitted or summation intensities, or
               an optimised combination of profile/sum."
      .expert_level = 1
    combine.Imid = None
      .type = floats
      .help = "A list of values to try for the midpoint, for profile/sum combination
               calculation: the value with the lowest Rmeas will be chosen.
               0 and 1 are special values that can be supplied to include profile
               and sum respectively in the comparison."
      .expert_level = 2
    combine.joint_analysis = True
      .type = bool
      .help = "Option of whether to do intensity combination optimisation
              separately (i.e. different Imid per dataset) or joint for
              multiple datasets"
      .expert_level = 2
  }
  weighting {
    weighting_scheme = *invvar unity GM cauchy huber
      .type = choice
      .help = "Weighting scheme used during Ih calculation. Weighting schemes
              other than invvar and unity may trigger iterative reweighting
              during minimisation, which may be unstable for certain minimisation
              engines (LBFGS)."
      .expert_level = 2
    optimise_errors = True
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
      .expert_level = 0
    error_model = *basic
      .type = choice
      .help = "The name of the error model to use, if optimise_errors is True."
      .expert_level = 1
    output_optimised_vars = True
      .type = bool
      .help = "If True, the error model determined will be applied to the
              intensity variances in the output files. This may result in
              a significant increase or decrease in the variances. The default
              is True as with the default inverse variance weighting scheme,
              the modified variances have been used as weights in scaling and
              therefore should be used as the variances when calculating merged
              intensities downstream. If this is distorting the data too much,
              then it is likely that the chosen error model is inappropriate."
      .expert_level = 2
  }
  cut_data {
    d_min = None
      .type = float
      .help = "Option to apply a high resolution cutoff for the dataset (i.e.
               the chosen reflections have d > d_min)."
      .expert_level = 1
    d_max = None
      .type = float
      .help = "Option to apply a low resolution cutoff for the dataset (i.e.
               the chosen reflections have d < d_max)."
      .expert_level = 1
    partiality_cutoff = 0.2
      .type = float
      .help = "Value below which reflections are removed from the dataset due
               to low partiality."
      .expert_level = 1
  }
  dataset_selection {
    use_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to use from a multi-dataset input."
      .expert_level = 2
    exclude_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to exclude from a multi-dataset input."
      .expert_level = 2
  }
  scaling_options {
    target_cycle = True
      .type = bool
      .help = "Option to turn of initial round of targeted scaling
               if some datasets are already scaled."
      .expert_level = 2
    only_target = False
      .type = bool
      .help = "Option to only do targeted scaling if some datasets
               are already scaled."
      .expert_level = 2
    only_save_targeted = True
      .type = bool
      .help = "If only_target is true, this option to change whether the dataset
              that is being scaled will be saved on its own, or combined with the
              already scaled dataset."
      .expert_level = 2
    target_model = None
      .type = path
      .help = "Path to cif file to use to calculate target intensities for
              scaling."
      .expert_level = 2
    target_mtz = None
      .type = path
      .help = "Path to merged mtz file to use as a target for scaling."
      .expert_level = 2
    nproc = 1
      .type = int(value_min=1)
      .help = "Number of blocks to divide the data into for minimisation.
              This also sets the number of processes to use if the option is
              available."
      .expert_level = 2
    use_free_set = False
      .type = bool
      .help = "Option to use a free set during scaling to check for overbiasing.
              This free set is used to calculate an RMSD, which is shown alongisde
              the 'working' RMSD during refinement, but is not currently used
              to terminate refinement or make any choices on the model."
      .expert_level = 2
    free_set_percentage = 10.0
      .type = float
      .help = "Percentage of symmetry equivalent groups to use for the free set,
              if use_free_set is True."
      .expert_level = 2
    free_set_offset = 0
      .type = int
      .help = "Offset for choosing unique groups for the free set from the whole
               set of unique groups."
      .expert_level = 2
    space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
      .expert_level = 1
    concurrent = True
      .type = bool
      .help = "Option to allow consecutive scaling if concurrent is
               set to False. The consecutive order is defined (and fixed)
               for each scaling model."
      .expert_level = 2
    full_matrix = True
      .type = bool
      .help = "Option to turn off GN/LM refinement round used to determine
               error estimates on scale factors."
      .expert_level = 2
    outlier_rejection = *standard simple
      .type = choice
      .help = "Choice of outlier rejection routine. Standard may take a
        significant amount of time to run for large datasets or high
        multiplicities, whereas simple should be quick for these datasets."
      .expert_level = 1
    outlier_zmax = 6.0
      .type = float(value_min=3.0)
      .help = "Cutoff z-score value for identifying outliers based on their
               normalised deviation within the group of equivalent reflections"
      .expert_level = 1
    verbosity = 2
      .type = int(value_min=0)
      .help = "The verbosity level"
      .expert_level = 2
  }

  ''')
