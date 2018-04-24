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
    E2_min = 0.8
      .type = float
      .help = "Minimum normalised E^2 value to select reflections for scaling"
    E2_max = 5.0
      .type = float
      .help = "Maximum normalised E^2 value to select reflections for scaling"
    Isigma_min = -5.0
      .type = float
      .help = "Option to use a I/sigma subset of reflections to determine scale factors"
    d_min = 0.0
      .type = float
      .help = "Option to use a d-value subset of reflections to determine scale factors"
  }
  weighting {
    weighting_scheme = *invvar tukey unity
      .type = choice
      .help = "Weighting scheme used during Ih calculation."
    optimise_error_model = False
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
    error_model_params = None
      .type = floats(size=2)
      .help = "Ability to force an error model adjustment, based on the model
              used in aimless - factors are called SDFac, SDadd in aimless."
  }
  cut_data {
    exclude_image_range = None
      .type = floats(size=2)
      .help = "Exclude a range of image numbers (start, stop) from the dataset,
               only used if a single dataset present."
    max_resolution = None
      .type = float
      .help = "Option to apply a maximum resolution cutoff for the dataset."
    min_resolution = None
      .type = float
      .help = "Option to apply a minimum resolution cutoff for the dataset."
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
    target_intensities = None
      .type = path
      .help = "Target intensity pickle file to use as target."
    space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
    concurrent = True
      .type = bool
      .help = "Option to allow consecutive scaling if concurrent is
               set to False"
    full_matrix_round = True
      .type = bool
      .help = "Option to turn off GN/LM refinement round used to determine
               error estimates on scale factors."
    outlier_rejection = standard *simple 0
      .type = choice
      .help = "Choice of outlier rejection routine. Standard may take a
        significant amount of time to run for large datasets or high
        multiplicities, whereas simple should be quick for these datasets."
    outlier_zmax = 9.0
      .type = float(value_min=6.0)
      .help = "Cutoff z-score value for identifying outliers based on their
               normalised deviation within the group of equivalent reflections"
    verbosity = 1
      .type = int(value_min=0)
      .help = "The verbosity level"
    integration_method = *prf sum
      .type = choice
      .help = "Option to choose from profile fitted or summation intensities."
  }

  ''')
