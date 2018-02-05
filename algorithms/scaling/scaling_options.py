import iotbx.phil

phil_scope = iotbx.phil.parse('''

  parameterisation {
    scale_term = True
      .type = bool
      .help = "Option to turn off decay correction (only for KB scaling)"
    scale_interval = 15.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the scale term"
    decay_term = True
      .type = bool
      .help = "Option to turn off decay correction"
    decay_interval = 20.0
      .type = float
      .help = "User specified rotation (phi) interval in degrees for phi binning
              for the decay term"
    absorption_term = True
      .type = bool
      .help = "Option to turn off absorption correction"
    lmax = 4
      .type = int
      .help = "Number of spherical harmonics to include for absorption correction,
              recommended to be no more than 6."
    surface_weight = 1e6
      .type = float
      .help = "Restraint weight applied to spherical harmonic terms in absorption
              correction."
    modulation_term = False
      .type = bool
      .help = "Option to turn off absorption correction"
  }
  reflection_selection {
    E2min = 0.8
      .type = float
      .help = "Minimum normalised E^2 value to select reflections for scaling"
    E2max = 5.0
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
    tukey_biweighting = False
      .type = bool
      .help = "Option to turn on tukey biweighting scheme for scaling weights."
    optimise_error_model = False
      .type = bool
      .help = "Option to allow optimisation of weights for scaling. Performs
               and additional scale factor minimisation after adjusting weights."
    error_model_params = None
      .type = floats(size=2)
      .help = "Ability to force an error model adjustment, using the model
              in aimless - factors are called SDFac, SDadd in aimless."
  }
  scaling_options {
    target = True
      .type = bool
      .help = "Option to turn of initial round of targeted scaling
               if some datasets are already scaled."
    only_target = False
      .type = bool
      .help = "Option to only do targeted scaling if some datasets
               are already scaled."
    force_space_group = None
      .type = str
      .help = "Option to specify space group for scaling"
    concurrent_scaling = True
      .type = bool
      .help = "Option to allow absorption correction after decay/scale, 
              if concurrent_scaling is set to False"
    reject_outliers = True
      .type = bool
      .help = "Option to turn on outlier rejection"
    outlier_zmax = 12.0
      .type = float
      .help = "Cutoff z-score value for identifying outliers based on their
               normalised deviation within the group of equivalent reflections"
    verbosity = 1
      .type = int(value_min=0)
      .help = "The verbosity level"
    integration_method = 'prf'
      .type = str
      .help = "Option to choose from profile fitted intensities (prf)
              or summation integrated intensities (sum)"
    minimisation_parameterisation = 'standard'
      .type = str
      .help = "Choice of 'standard' (multiplicative) or 'log' g-value
               minimisation parameterisation"
  }

  ''')
