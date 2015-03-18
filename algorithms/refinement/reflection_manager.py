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
from __future__ import division

from math import pi
from logging import info, debug

from scitbx import matrix
from dials.array_family import flex
from dials.algorithms.refinement import weighting_strategies

# constants
RAD2DEG = 180. / pi
DEG2RAD = pi / 180.

# helper functions
def calculate_entering_flags(reflections, experiments):
  """calculate entering flags for all reflections, and set them as a column
  of the reflection table."""

  # Init entering flags. These are always False for experiments that have no
  # rotation axis.
  enterings = flex.bool(len(reflections), False)

  for iexp, exp in enumerate(experiments):
    gonio = exp.goniometer
    if not gonio: continue
    axis = matrix.col(gonio.get_rotation_axis())
    s0 = matrix.col(exp.beam.get_s0())
    # calculate a unit vector normal to the spindle-beam plane for this
    # experiment, such that the vector placed at the centre of the Ewald sphere
    # points to the hemisphere in which reflections cross from inside to outside
    # of the sphere (reflections are exiting). NB this vector is in +ve Y
    # direction when using imgCIF coordinate frame.
    vec = s0.cross(axis)
    sel = reflections['id'] == iexp
    to_update = reflections['s1'].select(sel).dot(vec) < 0.
    enterings.set_selected(sel, to_update)

  return enterings

class BlockCalculator(object):
  """Utility class to calculate and set columns in the provided reflection
  table, which will be used during scan-varying refinement. The columns are a
  'block' number and an associated 'block_centre', giving the image number in
  the centre of the block"""

  def __init__(self, experiments, reflections):

    self._experiments = experiments
    self._reflections = reflections

    # do not create block column in the reflection table yet, in case we don't
    # need it at all

    return

  def _create_block_columns(self):
    """Create a column to contain the block number."""

    from scitbx.array_family import flex
    self._reflections['block'] = flex.size_t(len(self._reflections))
    self._reflections['block_centre'] = flex.double(len(self._reflections))
    return

  def per_width(self, width, deg=True):
    """Set blocks for all experiments according to a constant width"""

    if deg: width *= DEG2RAD
    self._create_block_columns()

    # get observed phi in radians
    phi_obs = self._reflections['xyzobs.mm.value'].parts()[2]

    for iexp, exp in enumerate(self._experiments):

      sel = self._reflections['id'] == iexp
      isel = sel.iselection()
      exp_phi = phi_obs.select(isel)

      start, stop = exp.scan.get_oscillation_range(deg=False)
      nblocks = int(abs(stop - start) / width) + 1
      # ensure width has the right sign and is wide enough that all reflections
      # get assigned a block
      _width = cmp(stop, start) * width + 1e-11
      half_width = width * (0.5 - 1e-11) # ensure round down behaviour

      block_starts = [start + n * _width for n in xrange(nblocks)]
      block_centres = [exp.scan.get_array_index_from_angle(
        e + half_width, deg=False) for e in block_starts]

      for b_num, (b_start, b_cent) in enumerate(zip(block_starts, block_centres)):
        sub_isel = isel.select((b_start <= exp_phi) & \
                                          (exp_phi <= (b_start + _width)))
        self._reflections['block'].set_selected(sub_isel, b_num)
        self._reflections['block_centre'].set_selected(sub_isel, b_cent)


    return self._reflections

  def per_image(self):
    """Set one block per image for all experiments"""

    self._create_block_columns()

    # get observed phi in radians
    phi_obs = self._reflections['xyzobs.mm.value'].parts()[2]

    for iexp, exp in enumerate(self._experiments):

      sel = self._reflections['id'] == iexp
      isel = sel.iselection()
      exp_phi = phi_obs.select(isel)

      # convert phi to integer frames
      frames = exp.scan.get_array_index_from_angle(exp_phi, deg=False)
      frames = flex.floor(frames).iround()

      start, stop = flex.min(frames), flex.max(frames)
      frame_range = range(start, stop + 1)

      for f_num, f in enumerate(frame_range):
        sub_isel = isel.select(frames == f)
        self._reflections['block'].set_selected(sub_isel, f_num)
        self._reflections['block_centre'].set_selected(sub_isel, f_num)

    return self._reflections


class ReflectionManager(object):
  """A class to maintain information about observed and predicted
  reflections for refinement.

  This new version keeps the reflections as a reflection table. Initialisation
  is not complete until the ReflectionManager is paired with a target function.
  Then, prediction can be done, followed by outlier rejection and any random
  sampling to form the working subset."""

  _weighting_strategy = weighting_strategies.StatisticalWeightingStrategy()

  def __init__(self, reflections,
                     experiments,
                     nref_per_degree=None,
                     min_num_obs=20,
                     max_sample_size=None,
                     min_sample_size=0,
                     close_to_spindle_cutoff=0.1,
                     iqr_multiplier=1.5,
                     weighting_strategy_override=None,
                     verbosity=0):

    # set verbosity
    # FIXME verbosity is deprecated now we use logging
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
    refs_to_keep = self._id_refs_to_keep(reflections)
    self._accepted_refs_size = len(refs_to_keep)

    # put full list of indexed reflections aside and select only the accepted
    # reflections to manage
    self._indexed = reflections
    self._reflections = reflections.select(flex.size_t(refs_to_keep))

    # keep minimum number of observations per experiment to allow as working set
    self._min_num_obs = min_num_obs

    # fail if there are too few reflections in the manager (the check also needs
    # to be repeated after outlier rejection and subsetting)
    self._check_too_few()

    # set entering flags for all reflections
    self._reflections['entering'] = calculate_entering_flags(self._reflections,
      self._experiments)

    # set observed frame numbers for all reflections if not already present
    self._calculate_frame_numbers()

    # set weights for all reflections
    if weighting_strategy_override is not None:
      self._weighting_strategy = weighting_strategy_override
    self._weighting_strategy.calculate_weights(self._reflections)

    # reset all use flags
    self.reset_accepted_reflections()

    # not known until the manager is finalised
    self._sample_size = None

    return

  def finalise(self):
    """Complete initialisation by performing outlier rejection and any
    requested subsetting. This function to be called by a Target object"""

    debug("Finalising the Reflection Manager")

    # print summary before outlier rejection
    self.print_stats_on_matches()

    # flag potential outliers
    rejection_occurred = self._flag_outliers()

    # delete all reflections from the manager that do not have a prediction
    # or were flagged as outliers
    msg = "Removing reflections not matched to predictions"
    if rejection_occurred: msg += " or marked as outliers"
    debug(msg)

    has_pred = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    inlier = ~self._reflections.get_flags(self._reflections.flags.centroid_outlier)
    self._reflections = self._reflections.select(has_pred & inlier)
    self._check_too_few()

    debug("%d reflections remain in the manager", len(self._reflections))

    # print summary after outlier rejection
    if rejection_occurred: self.print_stats_on_matches()

    # form working and free subsets
    self._create_working_set()
    self._sample_size = len(self._reflections)

    self._check_too_few()

    debug("Working set size = %d observations", self.get_sample_size())

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

  def _calculate_frame_numbers(self):
    """calculate observed frame numbers for all reflections, if not already
    set"""

    if self._reflections.has_key('xyzobs.px.value'): return

    # Frames are not set, so set them, with dummy observed pixel values
    frames = flex.double(len(self._reflections), 0.)
    for iexp, exp in enumerate(self._experiments):
      scan = exp.scan
      if not scan: continue
      sel = self._reflections['id'] == iexp
      xyzobs = self._reflections["xyzobs.mm.value"].select(sel)
      angles = xyzobs.parts()[2]
      to_update = scan.get_array_index_from_angle(angles, deg=False)
      frames.set_selected(sel, to_update)
    self._reflections['xyzobs.px.value'] = flex.vec3_double(
            flex.double(len(self._reflections), 0.),
            flex.double(len(self._reflections), 0.),
            frames)

    return

  def _id_refs_to_keep(self, obs_data):
    """Create a selection of observations that pass certain conditions.

    This step includes rejection of reflections too close to the spindle, and
    rejection of the (0,0,0) Miller index. Outlier rejection is done later."""

    # first exclude reflections with miller index set to 0,0,0
    sel = obs_data['miller_index'] != (0,0,0)
    inc = flex.size_t_range(len(obs_data)).select(sel)
    obs_data = obs_data.select(sel)

    # Default to True to pass test if there is no rotation axis for a particular
    # experiment
    to_keep = flex.bool(len(inc), True)

    for iexp, exp in enumerate(self._experiments):
      axis = self._axes[iexp]
      if not axis: continue
      sel = obs_data['id'] == iexp
      s0 = self._s0vecs[iexp]
      s1 = obs_data['s1'].select(sel)
      to_update = self._inclusion_test(s1, axis, s0)
      to_keep.set_selected(sel, to_update)

    inc = inc.select(to_keep)

    return inc

  def _inclusion_test(self, s1, axis, s0):
    """Test scattering vector s1 for inclusion"""

    # reject reflections for which the parallelepiped formed between
    # the gonio axis, s0 and s1 has a volume of less than the cutoff.
    # Those reflections are by definition closer to the spindle-beam
    # plane and for low values of the cutoff are troublesome to
    # integrate anyway.

    p_vol = flex.abs(s1.cross(flex.vec3_double(s1.size(), s0)).dot(axis))
    test = p_vol > self._close_to_spindle_cutoff

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
                    sweep_range_rad[0]) * RAD2DEG
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

    # create subsets
    free_sel = flex.bool(len(self._reflections), True)
    free_sel.set_selected(working_isel, False)
    self._free_reflections = self._reflections.select(free_sel)
    self._reflections = self._reflections.select(working_isel)

    return

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
    import copy
    sort_obs = copy.deepcopy(obs)
    if angular:
      sort_obs.sort('phi_resid', reverse=True)
    else:
      sort_obs['key'] = sort_obs['x_resid']**2 + sort_obs['y_resid']**2
      sort_obs.sort('key', reverse=True)
      del sort_obs['key']
    return sort_obs

  def get_indexed(self):
    """Return the reflections passed in as input"""

    return self._indexed

  def get_matches(self):
    """For every observation used in refinement return (a copy of) all data"""

    return self._reflections.select(self._reflections.get_flags(
      self._reflections.flags.used_in_refinement))

  def get_free_reflections(self):
    """Return all reflections that were accepted for refinement but not chosen
    in the working set"""

    return self._free_reflections

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    from libtbx.table_utils import simple_table
    from scitbx.math import five_number_summary
    x_resid = l['x_resid']
    y_resid = l['y_resid']
    phi_resid = l['phi_resid']
    w_x, w_y, w_phi = l['xyzobs.mm.weights'].parts()

    info("\nSummary statistics for observations matched to predictions:")
    header = ["", "Min", "Q1", "Med", "Q3", "Max"]
    rows = []
    row_data = five_number_summary(x_resid)
    rows.append(["Xc - Xo (mm)"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(y_resid)
    rows.append(["Yc - Yo (mm)"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(phi_resid)
    rows.append(["Phic - Phio (deg)"] + ["%.4g" % (e * RAD2DEG) for e in row_data])
    row_data = five_number_summary(w_x)
    rows.append(["X weights"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(w_y)
    rows.append(["Y weights"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(w_phi)
    rows.append(["Phi weights"] + ["%.4g" % (e * DEG2RAD**2) for e in row_data])
    st = simple_table(rows, header)
    info(st.format())
    info("")

    if len(l) < 20:
      debug("Fewer than 20 reflections matched!")
      return

    sl = self._sort_obs_by_residual(l)
    debug("Reflections with the worst 20 positional residuals:")
    debug("H, K, L, x_resid, y_resid, phi_resid, panel, x_obs, y_obs, " + \
      "phi_obs, weight_x_obs, weight_y_obs, weight_phi_obs")
    fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %d %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
    for i in xrange(20):
      e = sl[i]
      x_obs, y_obs, phi_obs = e['xyzobs.mm.value']
      msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                       e['y_resid'],
                       e['phi_resid'] * RAD2DEG,
                       e['panel'],
                       x_obs,
                       y_obs,
                       phi_obs * RAD2DEG,
                       e['xyzobs.mm.weights'][0],
                       e['xyzobs.mm.weights'][1],
                       e['xyzobs.mm.weights'][2] * DEG2RAD**2))
      debug(msg)
    sl = self._sort_obs_by_residual(sl, angular=True)
    debug("\nReflections with the worst 20 angular residuals:")
    debug("H, K, L, x_resid, y_resid, phi_resid, panel, x_obs, y_obs, " + \
      "phi_obs, weight_x_obs, weight_y_obs, weight_phi_obs")
    fmt = "(%3d, %3d, %3d) %5.3f %5.3f %6.4f %d %5.3f %5.3f %6.4f %5.3f %5.3f %6.4f"
    for i in xrange(20):
      e = sl[i]
      x_obs, y_obs, phi_obs = e['xyzobs.mm.value']
      msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                                             e['y_resid'],
                                             e['phi_resid'] * RAD2DEG,
                                             e['panel'],
                                             x_obs,
                                             y_obs,
                                             phi_obs * RAD2DEG,
                                             e['xyzobs.mm.weights'][0],
                                             e['xyzobs.mm.weights'][1],
                                             e['xyzobs.mm.weights'][2] * DEG2RAD**2))
      debug(msg)
    debug("")

    return

  def _flag_outliers(self):
    """Set the centroid outlier flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    sel = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    matches = self._reflections.select(sel)
    all_imatches = sel.iselection()

    max_panel = flex.max(matches['panel'])

    nreject = 0
    for pnl in range(max_panel + 1):
      sub_sel = matches['panel'] == pnl
      sub_matches = matches.select(sub_sel)
      sub_imatches = all_imatches.select(sub_sel)

      if len(sub_matches) >= self._min_num_obs:
        x_resid = sub_matches['x_resid']
        y_resid = sub_matches['y_resid']
        phi_resid = sub_matches['phi_resid']

        min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
        min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)
        min_phi, q1_phi, med_phi, q3_phi, max_phi = five_number_summary(phi_resid)

        iqr_x = q3_x - q1_x
        iqr_y = q3_y - q1_y
        iqr_phi = q3_phi - q1_phi

        cut_x = self._iqr_multiplier * iqr_x
        cut_y = self._iqr_multiplier * iqr_y
        cut_phi = self._iqr_multiplier * iqr_phi

        # accumulate a selection of outliers for this panel
        sel_out = x_resid > q3_x + cut_x
        sel_out.set_selected(x_resid < q1_x - cut_x, True)
        sel_out.set_selected(y_resid > q3_y + cut_y, True)
        sel_out.set_selected(y_resid < q1_y - cut_y, True)
        sel_out.set_selected(phi_resid > q3_phi + cut_phi, True)
        sel_out.set_selected(phi_resid < q1_phi - cut_phi, True)

        # get positions of outliers from the original matches
        ioutliers = sub_imatches.select(sel_out)
      else:
        debug("Only %d reflections on panel %d. All of these flagged " + \
              "as possible outliers.",
          len(sub_matches), pnl)
        ioutliers = sub_imatches

      # set those reflections as outliers
      self._reflections.set_flags(ioutliers,
        self._reflections.flags.centroid_outlier)

      nreject += len(ioutliers)

    if nreject == 0: return False

    info("%d reflections have been flagged as outliers", nreject)

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


class StillsReflectionManager(ReflectionManager):
  """Overloads for a Reflection Manager that does not exclude
  reflections too close to the spindle, and reports only information
  about X, Y, DelPsi residuals"""

  _weighting_strategy = weighting_strategies.StillsWeightingStrategy()

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()

    from libtbx.table_utils import simple_table
    from scitbx.math import five_number_summary
    x_resid = l['x_resid']
    y_resid = l['y_resid']
    delpsi = l['delpsical.rad']
    w_x, w_y, _ = l['xyzobs.mm.weights'].parts()
    w_delpsi = l['delpsical.weights']

    info("Summary statistics for observations matched to predictions:")
    header = ["", "Min", "Q1", "Med", "Q3", "Max"]
    rows = []
    row_data = five_number_summary(x_resid)
    rows.append(["Xc - Xo (mm)"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(y_resid)
    rows.append(["Yc - Yo (mm)"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(delpsi)
    rows.append(["DeltaPsi (deg)"] + ["%.4g" % (e * RAD2DEG) for e in row_data])
    row_data = five_number_summary(w_x)
    rows.append(["X weights"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(w_y)
    rows.append(["Y weights"] + ["%.4g" % e for e in row_data])
    row_data = five_number_summary(w_delpsi)
    rows.append(["DeltaPsi weights"] + ["%.4g" % (e * DEG2RAD**2) for e in row_data])
    st = simple_table(rows, header)
    info(st.format())

    if len(l) < 20:
      debug("Fewer than 20 reflections matched!")
      return

    sl = self._sort_obs_by_residual(l)
    debug("Reflections with the worst 20 positional residuals:")
    debug("H, K, L, x_resid, y_resid, weight_x_obs, weight_y_obs")
    fmt = "(%3d, %3d, %3d) %5.3f %5.3f %5.3f %5.3f"
    for i in xrange(20):
      e = sl[i]
      msg = fmt % tuple(e['miller_index'] + (e['x_resid'],
                                             e['y_resid'],
                                             e['xyzobs.mm.weights'][0],
                                             e['xyzobs.mm.weights'][1]))
      debug(msg)
    debug("\n")

    return

  def _flag_outliers(self):
    """Set the centroid outlier flag on matches with extreme (outlying) residuals.

    Outlier detection finds values more than _iqr_multiplier times the
    interquartile range from the quartiles. When x=1.5, this is Tukey's rule.

    Return boolean whether rejection was performed or not"""

    # return early if outlier rejection is disabled
    if self._iqr_multiplier is None: return False

    from scitbx.math import five_number_summary
    sel = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    matches = self._reflections.select(sel)
    all_imatches = sel.iselection()

    max_panel = flex.max(matches['panel'])

    nreject = 0
    for pnl in range(max_panel + 1):
      sub_sel = matches['panel'] == pnl
      sub_matches = matches.select(sub_sel)
      sub_imatches = all_imatches.select(sub_sel)

      if len(sub_matches) >= self._min_num_obs:
        x_resid = sub_matches['x_resid']
        y_resid = sub_matches['y_resid']
        delpsi = sub_matches['delpsical.rad']

        min_x, q1_x, med_x, q3_x, max_x = five_number_summary(x_resid)
        min_y, q1_y, med_y, q3_y, max_y = five_number_summary(y_resid)
        min_p, q1_p, med_p, q3_p, max_p = five_number_summary(delpsi)

        iqr_x = q3_x - q1_x
        iqr_y = q3_y - q1_y
        iqr_p = q3_p - q1_p

        cut_x = self._iqr_multiplier * iqr_x
        cut_y = self._iqr_multiplier * iqr_y
        cut_p = self._iqr_multiplier * iqr_p

        # accumulate a selection of outliers for this panel
        sel_out = x_resid > q3_x + cut_x
        sel_out.set_selected(x_resid < q1_x - cut_x, True)
        sel_out.set_selected(y_resid > q3_y + cut_y, True)
        sel_out.set_selected(y_resid < q1_y - cut_y, True)
        sel_out.set_selected(delpsi > q3_p + cut_p, True)
        sel_out.set_selected(delpsi < q1_p - cut_p, True)

        # get positions of outliers from the original matches
        ioutliers = sub_imatches.select(sel_out)
      else:
        debug("Only %d reflections on panel %d. All of these flagged " + \
              "as possible outliers.",
          len(sub_matches), pnl)
        ioutliers = sub_imatches

      # set those reflections as outliers
      self._reflections.set_flags(ioutliers,
        self._reflections.flags.centroid_outlier)

      nreject += len(ioutliers)

    if nreject == 0: return False

    info("%d reflections have been flagged as outliers" % nreject)

    return True

