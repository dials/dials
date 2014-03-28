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

from math import pi

from scitbx import matrix
from dials.array_family import flex

# constants
RAD_TO_DEG = 180. / pi

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
                     max_sample_size=None,
                     min_sample_size=0,
                     close_to_spindle_cutoff=0.1,
                     iqr_multiplier=1.5,
                     verbosity=0):

    # set verbosity
    self._verbosity = verbosity

    # keep track of models
    self._experiments = experiments
    goniometers = [e.goniometer for e in self._experiments]
    self._axes = [matrix.col(g.get_rotation_axis()) if g else None for g in goniometers]
    self._s0vecs = [matrix.col(e.beam.get_s0()) for e in self._experiments]

    # keep track of the original indices of the reflections
    reflections['iobs'] = flex.size_t_range(len(reflections))

    # set up the reflection inclusion cutoffs
    self._close_to_spindle_cutoff = close_to_spindle_cutoff #too close to spindle
    self._iqr_multiplier = iqr_multiplier #outlier rejection
    self._nref_per_degree = nref_per_degree #random subsets
    self._max_sample_size = max_sample_size #sample size ceiling
    self._min_sample_size = min_sample_size #sample size floor

    # exclude reflections that fail some inclusion criteria (currently just
    # close to spindle)
    self._input_size = len(reflections)
    refs_to_keep = self._id_refs_to_keep(reflections)
    self._accepted_refs_size = len(refs_to_keep)

    # select only the accepted reflections to manage
    self._reflections = reflections.select(flex.size_t(refs_to_keep))

    # keep minimum number of observations to allow as working set
    self._min_num_obs = min_num_obs

    # fail if there are too few reflections in the manager (the check also needs
    # to be repeated after outlier rejection and subsetting)
    self._check_too_few()

    # set entering flags for all reflections
    self._calculate_entering_flags()

    # set observed frame numbers for all reflections
    self._calculate_frame_numbers()

    # set weights for all reflections
    self._calculate_weights()

    # reset all use flags
    self.reset_accepted_reflections()

    # not known until the manager is finalised
    self._sample_size = None

    return

  def finalise(self):
    """Complete initialisation by performing outlier rejection and any
    requested subsetting. This function to be called by a Target object"""

    if self._verbosity > 1: print "Finalising the Reflection Manager"

    # print summary before outlier rejection
    self.print_stats_on_matches()

    # flag potential outliers
    rejection_occurred = self._reject_outliers()

    # delete all reflections from the manager that do not have a prediction
    # or were flagged as outliers
    if self._verbosity > 1:
      msg = "Removing reflections not matched to predictions"
      if rejection_occurred: msg += " or marked as outliers"
      print msg

    self._reflections = self._reflections.select(self._reflections.get_flags(
                      self._reflections.flags.used_in_refinement))
    self._check_too_few()

    if self._verbosity > 1:
      print len(self._reflections), "reflections remain in the manager"

    # print summary after outlier rejection
    if rejection_occurred: self.print_stats_on_matches()

    # form working subset
    refs_to_keep = self._create_working_set()
    self._sample_size = len(self._reflections)

    self._check_too_few()

    if self._verbosity > 1:
      print "Working set size = %d observations" % self.get_sample_size()

    return

  def _check_too_few(self):
    # fail if any of the experiments has too few reflections
    for iexp in range(len(self._experiments)):
      nref = (self._reflections['id'] == iexp).count(True)
      if nref < self._min_num_obs:
        msg = ('Remaining number of reflections = {0}, for experiment {1}, ' + \
               'which is below the configured limit for this reflection ' + \
               'manager').format(nref, iexp)
        raise RuntimeError(msg)
    return

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

    # Set entering flags. These are always False for experiments that have no
    # rotation axis.
    enterings = [matrix.col(ref['s1']).dot(vecs[ref['id']]) < 0. \
                 if vecs[ref['id']] else False for ref in self._reflections]

    self._reflections['entering'] = flex.bool(enterings)

    return

  def _calculate_frame_numbers(self):
    """calculate observed frame numbers for all reflections, if not already
    set"""

    scans = [e.scan for e in self._experiments]
    frames = [scans[ref['id']].get_array_index_from_angle(
                ref["xyzobs.mm.value"][2], deg=False) \
              if scans[ref['id']] else 0.0 for ref in self._reflections]
    #TEMP Test: if frames are already set, see that we have the same values
    # FIXME Remove the test if is is seen to pass in practice, and just return
    # early
    if self._reflections.has_key('xyzobs.px.value'):
      from libtbx.test_utils import approx_equal
      for a, b in zip(frames, self._reflections['xyzobs.px.value']):
        assert approx_equal(a, b[2])
    else: # Frames are not set, so set them, with dummy observed pixel values
      xyzobs = [(0., 0., f) for f in frames]
      self._reflections['xyzobs.px.value'] = flex.vec3_double(xyzobs)

    return

  def _calculate_weights(self):
    """set 'statistical weights', that is w(x) = 1/var(x)"""
    # TODO Plug-in to allow different weighting schemes

    weights = (self._reflections['xyzobs.mm.variance']).deep_copy()
    parts = weights.parts()
    for w in parts:
      sel = w >= 0.
      w.set_selected(sel, 1./w)
    self._reflections['xyzobs.mm.weights'] = flex.vec3_double(*parts)

    return

  def _id_refs_to_keep(self, obs_data):
    """Create a selection of observations that pass certain conditions.

    This step includes rejection of reflections too close to the spindle, and
    rejection of the (0,0,0) Miller index. Outlier rejection is done later."""

    #TODO Should be possible to 'vectorise' this using flex array operations
    #FIXME Allow inclusion test to pass if there is no rotation axis for
    # a particular experiment

    inc = [i for i, ref in enumerate(obs_data) if ref['miller_index'] != (0,0,0) \
           and self._inclusion_test(matrix.col(ref['s1']),
                                    self._axes[ref['id']],
                                    self._s0vecs[ref['id']])]

    return inc

  def _inclusion_test(self, s1, axis, s0):
    """Test scattering vector s1 for inclusion"""

    # reject reflections for which the parallelepiped formed between
    # the gonio axis, s0 and s1 has a volume of less than the cutoff.
    # Those reflections are by definition closer to the spindle-beam
    # plane and for low values of the cutoff are troublesome to
    # integrate anyway.

    test = abs(axis.dot(matrix.col(s1).cross(s0))) > \
        self._close_to_spindle_cutoff

    return test

  def _create_working_set(self):
    """Make a subset of the indices of reflections to use in refinement"""

    working_isel = flex.size_t()
    for iexp, exp in enumerate(self._experiments):

      sel = self._reflections['id'] == iexp
      isel = sel.iselection()
      #refs = self._reflections.select(sel)
      nrefs = sample_size = len(isel)

      # set sample size according to nref_per_degree (per experiment)
      if exp.scan and self._nref_per_degree:
        sweep_range_rad = exp.scan.get_oscillation_range(deg=False)
        width = abs(sweep_range_rad[1] -
                    sweep_range_rad[0]) * RAD_TO_DEG
        sample_size = int(self._nref_per_degree * width)
      else: sweep_range_rad = None

      # adjust sample size if below the chosen limit
      sample_size = max(sample_size, self._min_sample_size)

      # set maximum sample size if requested
      if self._max_sample_size:
        sample_size = min(sample_size, self._max_sample_size)

      # determine subset and collect indices
      if sample_size < nrefs:
        isel = isel.select(flex.random_selection(nrefs, sample_size))
      working_isel.extend(isel)

    # create subset
    self._reflections = self._reflections.select(working_isel)

    return

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

    sort_obs = deepcopy(obs)
    if angular:
      sort_obs.sort('phi_resid', reverse=True)
    else:
      sort_obs['key'] = sort_obs['x_resid']**2 + sort_obs['y_resid']**2
      sort_obs.sort('key', reverse=True)
      del sort_obs['key']
    return sort_obs

  def get_matches(self):
    """For every observation used in refinement return (a copy of) all data"""

    return self._reflections.select(self._reflections.get_flags(
      self._reflections.flags.used_in_refinement))

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    if self._verbosity > 1:

      from scitbx.math import five_number_summary
      x_resid = l['x_resid']
      y_resid = l['y_resid']
      phi_resid = l['phi_resid']
      w_x, w_y, w_phi = l['xyzobs.mm.weights'].parts()

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
          msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                           e['y_resid'],
                           e['phi_resid'],
                           e['weight_x_obs'],
                           e['weight_y_obs'],
                           e['weight_phi_obs']))
          print msg
        print
        sl = self._sort_obs_by_residual(sl, angular=True)
        print "\nReflections with the worst 20 angular residuals:"
        print "H, K, L, x_resid, y_resid, phi_resid, weight_x_obs, weight_y_obs, " + \
              "weight_phi_obs"
        fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
        for i in xrange(20):
          e = sl[i]
          msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                                                 e['y_resid'],
                                                 e['phi_resid'],
                                                 e['weight_x_obs'],
                                                 e['weight_y_obs'],
                                                 e['weight_phi_obs']))
          print msg
        print

    return

  def _reject_outliers(self):
    """Unset the use flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    sel = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    matches = self._reflections.select(sel)
    imatches = sel.iselection()

    x_resid = matches['x_resid']
    y_resid = matches['y_resid']
    phi_resid = matches['phi_resid']

    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
    min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)
    min_phi, q1_phi, med_phi, q3_phi, max_phi = five_number_summary(phi_resid)

    iqr_x = q3_x - q1_x
    iqr_y = q3_y - q1_y
    iqr_phi = q3_phi - q1_phi

    cut_x = self._iqr_multiplier * iqr_x
    cut_y = self._iqr_multiplier * iqr_y
    cut_phi = self._iqr_multiplier * iqr_phi

    # accumulate a selection of outliers
    sel = matches['x_resid'] > q3_x + cut_x
    sel.set_selected((matches['x_resid'] < q1_x - cut_x), True)
    sel.set_selected((matches['y_resid'] > q3_y + cut_y), True)
    sel.set_selected((matches['y_resid'] < q1_y - cut_y), True)
    sel.set_selected((matches['phi_resid'] > q3_phi + cut_phi), True)
    sel.set_selected((matches['phi_resid'] < q1_phi - cut_phi), True)

    # get positions of outliers from the original matches
    ioutliers = imatches.select(sel)

    # set those reflections to not be used
    self._reflections.unset_flags(ioutliers,
      self._reflections.flags.used_in_refinement)

    nreject = sel.count(True)
    if nreject == 0: return False

    if self._verbosity > 1:
      print "%d reflections have been rejected as outliers" % nreject

    return True


  def reset_accepted_reflections(self):
    """Reset use flags for all observations in preparation for a new set of
    predictions"""

    mask = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    self._reflections.unset_flags(mask, self._reflections.flags.used_in_refinement)
    return

  def get_obs(self):
    """Get the list of managed observations"""

    return self._reflections


class ReflectionManagerXY(ReflectionManager):
  """Overloads for a Reflection Manager that does not exclude
  reflections too close to the spindle, and reports only information
  about X, Y residuals"""

  def _id_refs_to_keep(self, obs_data):
    """For this version of the class, only reject the (0,0,0) reflections.
    We don't want to exclude reflections close to the spindle, as the spindle
    may not exist"""

    #FIXME Should not need to overload, if original version is fixed to
    # return true for inclusion whenever there is no rotation axis
    inc = [i for i, ref in enumerate(obs_data) if ref['miller_index'] != (0,0,0)]

    return inc

  # No need to overload the following. The nref_per_degree sampling won't be
  # done anyway if there is no scan
  # def _create_working_set(self):

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    if self._verbosity > 1:

      from scitbx.math import five_number_summary
      x_resid = l['x_resid']
      y_resid = l['y_resid']
      w_x, w_y, _ = l['xyzobs.mm.weights'].parts()

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
          msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                                                 e['y_resid'],
                                                 e['weight_x_obs'],
                                                 e['weight_y_obs']))
          print msg
        print

    return

  def _reject_outliers(self):
    """Unset the use flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    matches = self._reflections.select(self._reflections.get_flags(
      self._reflections.flags.used_in_refinement))
    imatches = matches.iselection()

    x_resid = matches['x_resid']
    y_resid = matches['y_resid']

    min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
    min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)

    iqr_x = q3_x - q1_x
    iqr_y = q3_y - q1_y

    cut_x = self._iqr_multiplier * iqr_x
    cut_y = self._iqr_multiplier * iqr_y

    # accumulate a selection of outliers
    sel = matches.select(matches['x_resid'] > q3_x + cut_x)
    sel.set_selected(matches.select(matches['x_resid'] < q1_x - cut_x))
    sel.set_selected(matches.select(matches['y_resid'] > q3_y + cut_y))
    sel.set_selected(matches.select(matches['y_resid'] < q1_y - cut_y))

    # get positions of outliers from the original matches
    ioutliers = imatches.select(sel)
    #mask = flex.bool(len(self._reflections))
    #mask.fill(False) # probably don't need this
    #mask.set_selected(ioutliers, True)

    # set those reflections to not be used
    self._reflections.unset_flags(ioutliers,
      self._reflections.flags.used_in_refinement)

    nreject = sel.count(True)
    if nreject == 0: return False

    if self._verbosity > 1:
      print "%d reflections have been rejected as outliers" % nreject

    return True
