#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""Contains classes used to manage the reflections used during refinement,
principally ReflectionManager."""

class ObsPredMatch:
  """
  A bucket class containing data for a prediction that has been
  matched to an observation.

  This contains all the raw material needed to calculate the target function
  value and gradients
  """

  # initialise with an observation
  def __init__(self, iobs, crystal, hkl, entering, frame, panel,
                     x_obs, sigx_obs, weight_x_obs,
                     y_obs, sigy_obs, weight_y_obs,
                     phi_obs, sigphi_obs, weight_phi_obs):

    self.iobs = iobs
    self.crystal_id = crystal
    self.miller_index = hkl
    self.entering = entering
    self.frame_obs = frame
    self.panel = panel
    self.x_obs = x_obs
    self.sigx_obs = sigx_obs
    self.weight_x_obs = weight_x_obs

    self.y_obs = y_obs
    self.sigy_obs = sigy_obs
    self.weight_y_obs = weight_y_obs

    self.phi_obs = phi_obs
    self.sigphi_obs = sigphi_obs
    self.weight_phi_obs = weight_phi_obs

    self.x_calc = None
    self.y_calc = None
    self.phi_calc = None
    self.s_calc = None

    # gradients will be a list, of length equal to the number of free
    # parameters, whose elements are the gradients in the space of
    # the residuals, e.g. (dX/dp, dY/dp, dPhi/dp)
    self.gradients = None

    self.y_resid = None
    self.y_resid2 = None
    self.x_resid = None
    self.x_resid2 = None
    self.phi_resid = None
    self.phi_resid2 = None

    self.is_matched = False

  # update with a prediction
  def update_prediction(self, x_calc, y_calc, phi_calc, s_calc, gradients):

    self.x_calc = x_calc
    self.y_calc = y_calc
    self.phi_calc = phi_calc
    self.s_calc = s_calc

    self.gradients = gradients

    # calculate residuals
    self.x_resid = x_calc - self.x_obs
    self.x_resid2 = self.x_resid**2
    self.y_resid = y_calc - self.y_obs
    self.y_resid2 = self.y_resid**2
    self.phi_resid = phi_calc - self.phi_obs
    self.phi_resid2 = self.phi_resid**2

    self.is_matched = True

  def reset(self):

    """Flag this observation to not be used"""
    self.is_matched = False


class ReflectionManager(object):
  """A class to maintain information about observed and predicted
  reflections for refinement.

  This new version keeps the reflections as a reflection table. Initialisation
  is not complete until the ReflectionManager is paired with a target function.
  Then, prediction can be done, followed by outlier rejection and any random
  sampling to form the working subset."""

  def __init__(self, reflections,
                     experiments,
                     nref_per_degree=None,
                     min_num_obs=20,
                     max_num_obs=None,
                     minimum_sample_size=0,
                     close_to_spindle_cutoff=0.1,
                     iqr_multiplier=1.5,
                     verbosity=0):

    # track whether this is the first update of predictions or not
    # DEPRECATED
    #self.first_update = True

    # set verbosity
    self._verbosity = verbosity

    ####FIXME currently extracting quantities from the FIRST Experiment only
    # DEPRECATED. KEEP THE EXPERIMENTS ENTIRE
    #if experiments[0].scan:
    #  sweep_range_rad = experiments[0].scan.get_oscillation_range(deg=False)
    #else: sweep_range_rad = None
    #gonio = experiments[0].goniometer
    #beam = experiments[0].beam
    ####
    self._experiments = experiments
    self._goniometers = [e.goniometer for e in self._experiments]
    self._axes = [matrix.col(g.get_rotation_axis()) if g else None for g in goniometers]
    self._s0vecs = [matrix.col(e.beam.get_s0()) for e in self.experiments]
    # keep references to the beam, goniometer and sweep range (for
    # reflection exclusion and subsetting)
    # DEPRECATED - USE THE EXPERIMENTS
    #self._beam = beam
    #self._gonio = gonio
    #self._sweep_range_rad = sweep_range_rad

    # find vector normal to the spindle-beam plane for the initial model
    #vecn = self._spindle_beam_plane_normal()

    # set up the reflection inclusion cutoffs
    self._close_to_spindle_cutoff = close_to_spindle_cutoff
    self._iqr_multiplier = iqr_multiplier

    # exclude reflections that fail inclusion criteria
    self._input_size = len(reflections)
    refs_to_keep = self._id_refs_to_keep(reflections)
    self._accepted_refs_size = len(refs_to_keep)

    # select only the accepted reflections to manage
    self._reflections = reflections.select(flex.size_t(refs_to_keep))

    # choose a random subset of data for refinement
    # DO THIS LATER, AFTER PREDICTION AND OUTLIER REJECTION
    #self._sample_size = self._accepted_refs_size
    #refs_to_keep = self._create_working_set(refs_to_keep,
    #                                        nref_per_degree,
    #                                        minimum_sample_size,
    #                                        max_num_obs)

    # keep minimum number of observations to allow as working set
    self._min_num_obs = min_num_obs

    # store observation information in a list of observation-prediction
    # pairs (prediction information will go in here later)
    #self._obs_pred_pairs = []

    # calculate entering flags for all reflections
    self._calculate_entering_flags()

    for i in refs_to_keep:

      ref = reflections[i]
      exp_id = ref['id']
      h = ref['miller_index']
      s = matrix.col(ref['s1'])
      entering = s.dot(vecn) < 0.
      panel = ref['panel']
      x = ref["xyzobs.mm.value"][0]
      y = ref["xyzobs.mm.value"][1]
      phi = ref["xyzobs.mm.value"][2]
      if experiments[exp_id].scan:
        frame = experiments[exp_id].scan.get_array_index_from_angle(phi, deg=False)
      else:
        frame = 0
      sig_x, sig_y, sig_phi = [sqrt(e) for e in ref["xyzobs.mm.variance"]]
      w_x = w_y = w_phi = 0
      if ref["xyzobs.mm.variance"][0] != 0: w_x   = 1./ref["xyzobs.mm.variance"][0]
      if ref["xyzobs.mm.variance"][1] != 0: w_y   = 1./ref["xyzobs.mm.variance"][1]
      if ref["xyzobs.mm.variance"][2] != 0: w_phi = 1./ref["xyzobs.mm.variance"][2]

      self._obs_pred_pairs.append(ObsPredMatch(i, exp_id, h,
                                               entering, frame, panel,
                                               x, sig_x, w_x,
                                               y, sig_y, w_y,
                                               phi, sig_phi, w_phi))

    # fail if there are too few reflections in the manager
    if len(self._obs_pred_pairs) < self._min_num_obs:
      msg = ('Remaining number of reflections = {0}, which is below '+ \
          'the configured limit for creating this reflection ' + \
          'manager').format(len(self._obs_pred_pairs))
      raise RuntimeError(msg)

  def _calculate_entering_flags(self):
    """calculate entering flags for all reflections, and set them as a column
    of the reflection table."""

    # calculate unit vectors normal to the spindle-beam plane for each
    # experiment, such that the vector placed at the centre of the Ewald sphere
    # points to the hemisphere in which reflections cross from inside to outside
    # of the sphere (reflections are exiting). NB this vector is in +ve Y
    # direction when using imgCIF coordinate frame.
    vecs = [self._s0vecs[i].cross(e).normalize() \
            if e else None for i, e in enumerate(self._axes)]

    enterings = [ref[''].dot(vecn) < 0.]

  def _id_refs_to_keep(self, obs_data):
    """Create a selection of observations that pass certain conditions.

    This step includes rejection of reflections too close to the spindle, and
    rejection of the (0,0,0) Miller index. Outlier rejection is done later."""



    inc = [i for i, ref in enumerate(obs_data) if ref['miller_index'] != (0,0,0) \
           and self._inclusion_test(matrix.col(ref['s1']),
                                    self._axes[ref['id']],
                                    self._s0vecs[ref['id']])]

    return inc

  def _inclusion_test(self, s1, axis, s0):
    """Test scattering vector s for inclusion"""

    # reject reflections for which the parallelepiped formed between
    # the gonio axis, s0 and s1 has a volume of less than the cutoff.
    # Those reflections are by definition closer to the spindle-beam
    # plane and for low values of the cutoff are troublesome to
    # integrate anyway.

    test = abs(axis.dot(matrix.col(s1).cross(s0))) > \
        self._close_to_spindle_cutoff

    return test

  def _create_working_set(self, indices, nref_per_degree,
                          minimum_sample_size, max_num_obs):
    """Make a subset of the indices of reflections to use in refinement"""

    working_indices = indices
    sample_size = len(working_indices)

    # set sample size according to nref_per_degree
    if nref_per_degree:
      width = abs(self._sweep_range_rad[1] -
                  self._sweep_range_rad[0]) * RAD_TO_DEG
      sample_size = int(nref_per_degree * width)

    # adjust if this is below the chosen limit
    sample_size = max(sample_size, minimum_sample_size)

    # set maximum sample size
    if max_num_obs:
      if sample_size > max_num_obs:
        sample_size = max_num_obs

    # sample the data and record the sample size
    if sample_size < len(working_indices):
      self._sample_size = sample_size
      working_indices = random.sample(working_indices,
                                      self._sample_size)
    return(working_indices)

  def get_input_size(self):
    """Return the number of observations in the initial list supplied
    as input"""

    return self._input_size

  def get_accepted_refs_size(self):
    """Return the number of observations that pass inclusion criteria and
    can potentially be used for refinement"""

    return self._accepted_refs_size

  def get_sample_size(self):
    """Return the number of observations in the working set to be
    used for refinement"""

    return self._sample_size

  def _sort_obs_by_residual(self, obs, angular=False):
    """For diagnostic purposes, sort the obs-pred matches so that the
    highest residuals are first. By default, sort by positional
    residual, unless angular=True.

    The earliest entries in the return list may be those that are
    causing problems in refinement.

    """

    if angular:
      k = lambda e: e.phi_resid
    else:
      k = lambda e: e.x_resid**2 + e.y_resid**2
    sort_obs = sorted(obs, key=k, reverse=True)

    return sort_obs

  def get_matches(self):
    """For every observation matched with a prediction return all data"""

    return [obs for obs in self._obs_pred_pairs if obs.is_matched]

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()
    if self._verbosity > 1:

      from scitbx.math import five_number_summary
      x_resid = [e.x_resid for e in l]
      y_resid = [e.y_resid for e in l]
      phi_resid = [e.phi_resid for e in l]
      w_x = [e.weight_x_obs for e in l]
      w_y = [e.weight_y_obs for e in l]
      w_phi = [e.weight_phi_obs for e in l]

      print "\nSummary statistics for observations matched to predictions:"
      print ("                      "
             "Min         Q1        Med         Q3        Max")
      print "(Xc-Xo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(x_resid))
      print "(Yc-Yo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(y_resid))
      print "(Phic-Phio)    {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(phi_resid))
      print "X weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_x))
      print "Y weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_y))
      print "Phi weights    {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_phi))
      print

      if len(l) >= 20 and self._verbosity > 2:

        sl = self._sort_obs_by_residual(l)
        print "Reflections with the worst 20 positional residuals:"
        print "H, K, L, x_resid, y_resid, phi_resid, weight_x_obs, weight_y_obs, " + \
              "weight_phi_obs"
        fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
        for i in xrange(20):
          e = sl[i]
          msg = fmt % tuple(e.miller_index + (e.x_resid,
                           e.y_resid,
                           e.phi_resid,
                           e.weight_x_obs,
                           e.weight_y_obs,
                           e.weight_phi_obs))
          print msg
        print
        sl = self._sort_obs_by_residual(sl, angular=True)
        print "\nReflections with the worst 20 angular residuals:"
        print "H, K, L, x_resid, y_resid, phi_resid, weight_x_obs, weight_y_obs, " + \
              "weight_phi_obs"
        fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
        for i in xrange(20):
          e = sl[i]
          msg = fmt % tuple(e.miller_index + (e.x_resid,
                           e.y_resid,
                           e.phi_resid,
                           e.weight_x_obs,
                           e.weight_y_obs,
                           e.weight_phi_obs))
          print msg
        print

    return

  def reject_outliers(self):
    """Unset the use flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    matches = [obs for obs in self._obs_pred_pairs if obs.is_matched]

    x_resid = [e.x_resid for e in matches]
    y_resid = [e.y_resid for e in matches]
    phi_resid = [e.phi_resid for e in matches]

    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
    min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)
    min_phi, q1_phi, med_phi, q3_phi, max_phi = five_number_summary(phi_resid)

    iqr_x = q3_x - q1_x
    iqr_y = q3_y - q1_y
    iqr_phi = q3_phi - q1_phi

    cut_x = self._iqr_multiplier * iqr_x
    cut_y = self._iqr_multiplier * iqr_y
    cut_phi = self._iqr_multiplier * iqr_phi

    for m in matches:
      if m.x_resid > q3_x + cut_x: m.is_matched = False
      if m.x_resid < q1_x - cut_x: m.is_matched = False
      if m.y_resid > q3_y + cut_y: m.is_matched = False
      if m.y_resid < q1_y - cut_y: m.is_matched = False
      if m.phi_resid > q3_phi + cut_phi: m.is_matched = False
      if m.phi_resid < q1_phi - cut_phi: m.is_matched = False

    nreject = [m.is_matched for m in matches].count(False)

    if nreject == 0: return False

    if self._verbosity > 1:
      print "%d reflections have been rejected as outliers" % nreject

    return True

  def strip_unmatched_observations(self):
    """Delete observations from the manager that are not matched to a
    prediction. Typically used once, after the first update of
    predictions."""

    if self._verbosity > 1:
      print "Removing reflections not matched to predictions"

    self._obs_pred_pairs = [e for e in self._obs_pred_pairs if e.is_matched]

    if len(self._obs_pred_pairs) < self._min_num_obs:
      msg = ('Remaining number of reflections = {0}, which is below '+ \
          'the configured limit for this reflection manager').format(
              len(self._obs_pred_pairs))
      raise RuntimeError(msg)

    if self._verbosity > 1:
      print len(self._obs_pred_pairs), "reflections remain in the manager"

    return

  def reset_accepted_reflections(self):
    """Reset all observations to use=False in preparation for a new set of
    predictions"""

    for obs in self._obs_pred_pairs: obs.reset()

  def get_obs(self):
    """Get the list of managed observations"""

    return self._obs_pred_pairs


class ReflectionManagerXY(ReflectionManager):
  """Overloads for a Reflection Manager that does not exclude
  reflections too close to the spindle, and reports only information
  about X, Y residuals"""

  def _spindle_beam_plane_normal(self):
    """There is no goniometer, so overload to return None"""

    return None

  def _id_refs_to_keep(self, obs_data):
    """For this version of the class, only reject the (0,0,0) reflections.
    We don't want to exclude reflections close to the spindle, as the spindle
    may not exist"""

    inc = [i for i, ref in enumerate(obs_data) if ref['miller_index'] != (0,0,0)]

    return inc

  def _create_working_set(self, indices, nref_per_degree,
                                         minimum_sample_size,
                                         max_num_obs):
    """Make a subset of the indices of reflections to use in refinement.

    This version ignores nref_per_degree"""

    working_indices = indices
    sample_size = len(working_indices)

    # set maximum sample size
    if max_num_obs:
      if sample_size > max_num_obs:
        sample_size = max_num_obs

    # sample the data and record the sample size
    if sample_size < len(working_indices):
      self._sample_size = sample_size
      working_indices = random.sample(working_indices,
                                      self._sample_size)
    return(working_indices)

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    if self._verbosity > 1:

      from scitbx.math import five_number_summary
      x_resid = [e.x_resid for e in l]
      y_resid = [e.y_resid for e in l]
      w_x = [e.weight_x_obs for e in l]
      w_y = [e.weight_y_obs for e in l]

      print "\nSummary statistics for observations matched to predictions:"
      print ("                      "
             "Min         Q1        Med         Q3        Max")
      print "(Xc-Xo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(x_resid))
      print "(Yc-Yo)        {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(y_resid))
      print "X weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_x))
      print "Y weights      {0:10.5g} {1:10.5g} {2:10.5g} {3:10.5g} {4:10.5g}".\
        format(*five_number_summary(w_y))
      print

      if len(l) >= 20 and self._verbosity > 2:

        sl = self._sort_obs_by_residual(l)
        print "Reflections with the worst 20 positional residuals:"
        print "H, K, L, x_resid, y_resid, weight_x_obs, weight_y_obs"
        fmt = "(%3d, %3d, %3d) %5.3f %5.3f %5.3f %5.3f"
        for i in xrange(20):
          e = sl[i]
          msg = fmt % tuple(e.miller_index + (e.x_resid,
                           e.y_resid,
                           e.weight_x_obs,
                           e.weight_y_obs))
          print msg
        print

    return

  def reject_outliers(self):
    """Unset the use flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    matches = [obs for obs in self._obs_pred_pairs if obs.is_matched]

    x_resid = [e.x_resid for e in matches]
    y_resid = [e.y_resid for e in matches]

    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
    min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)

    iqr_x = q3_x - q1_x
    iqr_y = q3_y - q1_y

    cut_x = self._iqr_multiplier * iqr_x
    cut_y = self._iqr_multiplier * iqr_y

    for m in matches:
      if m.x_resid > q3_x + cut_x: m.is_matched = False
      if m.x_resid < q1_x - cut_x: m.is_matched = False
      if m.y_resid > q3_y + cut_y: m.is_matched = False
      if m.y_resid < q1_y - cut_y: m.is_matched = False

    nreject = [m.is_matched for m in matches].count(False)

    if nreject == 0: return False

    if self._verbosity > 1:
      print "%d reflections have been rejected as outliers" % nreject

    return True
