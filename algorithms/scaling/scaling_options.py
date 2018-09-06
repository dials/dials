"""
Phil scope of options for scaling.
"""
import iotbx.phil

phil_scope = iotbx.phil.parse('''

  parameterisation {
    scale_term = True
      .type = bool
      .help = "Option to turn off scale correction (for physical/KB
               default models)."
    scale_interval = 15.0
      .type = float(value_min=1.0)
      .help = "Rotation (phi) interval between model parameters for the scale
               component (physical model)."
    decay_term = True
      .type = bool
      .help = "Option to turn off decay correction (for physical/array/KB
               default models)."
    decay_interval = 20.0
      .type = float(value_min=1.0)
      .help = "Rotation (phi) interval between model parameters for the decay
               component (physical/array default models)."
    n_resolution_bins = 10
      .type = int(value_min=1)
      .help = "Number of resolution bins to use for the decay term in the
               array-based model."
    absorption_term = True
      .type = bool
      .help = "Option to turn off absorption correction (for physical/array
               default models)."
    lmax = 4
      .type = int(value_min=2)
      .help = "Number of spherical harmonics to include for absorption
              correction (for physical default model), recommended to be no
              more than 6."
    surface_weight = 1e6
      .type = float(value_min=0.0)
      .help = "Restraint weight applied to spherical harmonic terms in the
               physical model absorption correction."
    modulation_term = False
      .type = bool
      .help = "Option to turn on a detector correction for the array default
               model."
    n_modulation_bins = 20
      .type = int(value_min=1)
      .help = "Number of bins in each dimension (applied to both x and y) for
              binning the detector position for the modulation term of the
              array model."
    n_absorption_bins = 3
      .type = int(value_min=1)
      .help = "Number of bins in each dimension (applied to both x and y) for
              binning the detector position for the absorption term of the
              array model."
  }
  reflection_selection {
    E2_range = 0.8, 5.0
      .type = floats(size=2)
      .help = "Minimum and maximum normalised E^2 value to used to select a
              subset of reflections for minimising the scaling model."
    Isigma_range = -5.0, 0.0
      .type = floats(size=2)
      .help = "Minimum and maximum I/sigma values used to subset of reflections
              to determine the scaling model. Here a value of 0.0 for the max
              means no limit applied."
    d_range = None
      .type = floats(size=2)
      .help = "Minimum and maximum - values used to subset of reflections
              to determine the scaling model."
    min_partiality = 0.6
      .type = float
      .help = "Minimum partiality to use when selecting subset of reflections
               to determine the scaling model."
    intensity_choice = profile sum *combine
      .type = choice
      .help = "Option to choose from profile fitted or summation intensities, or
               an optimised combination of profile/sum."
    combine.Imid = None
      .type = floats
      .help = "A list of values to try for the midpoint, for profile/sum combination
               calculation: the value with the lowest Rmeas will be chosen.
               0 and 1 are special values that can be supplied to include profile
               and sum respectively in the comparison."
    combine.joint_analysis = True
      .type = bool
      .help = "Option of whether to do intensity combination optimisation
              separately (i.e. different Imid per dataset) or joint for
              multiple datasets"
  }
  weighting {
    weighting_scheme = *invvar unity GM cauchy huber
      .type = choice
      .help = "Weighting scheme used during Ih calculation. Weighting schemes
              other than invvar and unity may trigger iterative reweighting
              during minimisation, which may be unstable for certain minimisation
              engines (LBFGS)."
    optimise_errors = True
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
    error_model = *basic
      .type = choice
      .help = "The name of the error model to use, if optimise_errors is True."
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
  }
  cut_data {
    exclude_image_range = None
      .type = floats(size=2)
      .help = "Exclude a range of image numbers (start, stop) from the dataset,
               only used if a single dataset present."
    d_min = None
      .type = float
      .help = "Option to apply a high resolution cutoff for the dataset (i.e.
               the chosen reflections have d > d_min)."
    d_max = None
      .type = float
      .help = "Option to apply a low resolution cutoff for the dataset (i.e.
               the chosen reflections have d < d_max)."
    partiality_cutoff = 0.2
      .type = float
      .help = "Value below which reflections are removed from the dataset due
               to low partiality."
  }
  dataset_selection {
    use_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to use from a multi-dataset input."
    exclude_datasets = None
      .type = strings
      .help = "Choose a subset of datasets, based on the dataset id (as defined
               in the reflection table), to exclude from a multi-dataset input."
  }
  scaling_options {
    target_cycle = True
      .type = bool
      .help = "Option to turn of initial round of targeted scaling
               if some datasets are already scaled."
    only_target = False
      .type = bool
      .help = "Option to only do targeted scaling if some datasets
               are already scaled."
    only_save_targeted = True
      .type = bool
      .help = "If only_target is true, this option to change whether the dataset
              that is being scaled will be saved on its own, or combined with the
              already scaled dataset."
    target_model = None
      .type = path
      .help = "Path to cif file to use to calculate target intensities for
              scaling."
    target_mtz = None
      .type = path
      .help = "Path to merged mtz file to use as a target for scaling."
    nproc = 1
      .type = int(value_min=1)
      .help = "Number of blocks to divide the data into for minimisation.
              This also sets the number of processes to use if the option is
              available."
    use_free_set = False
      .type = bool
      .help = "Option to use a free set during scaling to check for overbiasing.
              This free set is used to calculate an RMSD, which is shown alongisde
              the 'working' RMSD during refinement, but is not currently used
              to terminate refinement or make any choices on the model."
    free_set_percentage = 10.0
      .type = float
      .help = "Percentage of symmetry equivalent groups to use for the free set,
              if use_free_set is True."
    free_set_offset = 0
      .type = int
      .help = "Offset for choosing unique groups for the free set from the whole
               set of unique groups."
    space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
    concurrent = True
      .type = bool
      .help = "Option to allow consecutive scaling if concurrent is
               set to False. The consecutive order is defined (and fixed)
               for each scaling model."
    full_matrix = True
      .type = bool
      .help = "Option to turn off GN/LM refinement round used to determine
               error estimates on scale factors."
    outlier_rejection = *standard simple
      .type = choice
      .help = "Choice of outlier rejection routine. Standard may take a
        significant amount of time to run for large datasets or high
        multiplicities, whereas simple should be quick for these datasets."
    outlier_zmax = 6.0
      .type = float(value_min=3.0)
      .help = "Cutoff z-score value for identifying outliers based on their
               normalised deviation within the group of equivalent reflections"
    verbosity = 2
      .type = int(value_min=0)
      .help = "The verbosity level"
  }

  ''')
