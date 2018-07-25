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
from __future__ import absolute_import, division

from math import pi
import logging
logger = logging.getLogger(__name__)

import libtbx
from scitbx import matrix
from dials.array_family import flex
from dials.algorithms.refinement import weighting_strategies
from dials.algorithms.refinement.analysis.centroid_analysis import \
  CentroidAnalyser
from dials.algorithms.refinement.refinement_helpers import \
  calculate_frame_numbers, set_obs_s1

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
        f_cent = f + 0.5
        self._reflections['block'].set_selected(sub_isel, f_num)
        self._reflections['block_centre'].set_selected(sub_isel, f_cent)

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
                     max_sample_size=None,
                     min_sample_size=0,
                     close_to_spindle_cutoff=0.02,
                     trim_scan_edges=0.0,
                     outlier_detector=None,
                     weighting_strategy_override=None,
                     verbosity=0):

    # set verbosity
    if verbosity == 0:
      logger.disabled = True
    self._verbosity = verbosity

    # keep track of models
    self._experiments = experiments
    goniometers = [e.goniometer for e in self._experiments]
    self._axes = [matrix.col(g.get_rotation_axis()) if g else None for g in goniometers]
    self._s0vecs = [matrix.col(e.beam.get_s0()) for e in self._experiments]

    # unset the refinement flags (creates flags field if needed)
    reflections.unset_flags(flex.size_t_range(len(reflections)),
        flex.reflection_table.flags.used_in_refinement)

    # check that the observed beam vectors are stored: if not, compute them
    n_s1_set = set_obs_s1(reflections, experiments)
    if n_s1_set > 0:
      logger.debug("Set scattering vectors for %d reflections", n_s1_set)

    # keep track of the original indices of the reflections
    reflections['iobs'] = flex.size_t_range(len(reflections))

    #Check for monotonically increasing value range. If not, ref_table isn't sorted,
    # and proceed to sort by id and panel. This is required for the C++ extension
    # modules to allow for nlogn subselection of values used in refinement.
    l_id = reflections["id"]
    id0 = l_id[0]
    for ii in xrange(1,len(l_id)):
      if id0 <= l_id[ii]:
        id0 = l_id[ii]
      else:
        reflections.sort("id") #Ensuring the ref_table is sorted by id
        reflections.subsort("id","panel") #Ensuring that within each sorted id block, sorting is next performed by panel
        break

    # set up the reflection inclusion criteria
    self._close_to_spindle_cutoff = close_to_spindle_cutoff # close to spindle
    self._trim_scan_edges = DEG2RAD * trim_scan_edges # close to the scan edge
    self._outlier_detector = outlier_detector # for outlier rejection
    self._nref_per_degree = nref_per_degree # random subsets
    self._max_sample_size = max_sample_size # sample size ceiling
    self._min_sample_size = min_sample_size # sample size floor

    # exclude reflections that fail some inclusion criteria
    refs_to_keep = self._id_refs_to_keep(reflections)
    self._accepted_refs_size = len(refs_to_keep)

    # set entering flags for all reflections
    reflections['entering'] = calculate_entering_flags(reflections,
      self._experiments)

    # set observed frame numbers for all reflections if not already present
    calculate_frame_numbers(reflections, self._experiments)

    # reset all use flags
    self.reset_accepted_reflections(reflections)

    # put full list of indexed reflections aside and select only the reflections
    # that were not excluded to manage
    self._indexed = reflections
    self._reflections = reflections.select(flex.size_t(refs_to_keep))

    # set weights for all kept reflections
    if weighting_strategy_override is not None:
      self._weighting_strategy = weighting_strategy_override
    self._weighting_strategy.calculate_weights(self._reflections)

    # not known until the manager is finalised
    self._sample_size = None

    return

  def get_centroid_analyser(self, debug=False):
    """Create a CentroidAnalysis object for the current reflections"""

    return CentroidAnalyser(self._reflections, debug=debug)

  def finalise(self, analysis=None):
    """Complete initialisation by performing outlier rejection and any
    requested subsetting. If a list of results from a CentroidAnalysis
    object is provided, these may be used to determine outlier rejection
    block widths"""

    logger.debug("Finalising the Reflection Manager")

    # print summary before outlier rejection
    if self._verbosity > 1: self.print_stats_on_matches()

    # reset centroid_outlier flags in both the working reflections and the
    # original indexed reflections
    mask = self._reflections.get_flags(self._reflections.flags.centroid_outlier)
    self._reflections.unset_flags(mask, self._reflections.flags.centroid_outlier)
    mask = self._indexed.get_flags(self._indexed.flags.centroid_outlier)
    self._indexed.unset_flags(mask, self._indexed.flags.centroid_outlier)

    # outlier rejection if requested
    if self._outlier_detector is None:
      rejection_occurred = False
    else:
      if self._outlier_detector.get_block_width() is libtbx.Auto:
        if analysis is None:
          # without analysis available, set 18.0 degrees universally
          self._outlier_detector.set_block_width(18.0)
        else:
          # with analysis, choose the maximum of 18 degrees or the block size
          # for each experiment
          widths = [e.get('block_size') for e in analysis]
          widths = [max(e, 18.0) if e is not None else None for e in widths]
          self._outlier_detector.set_block_width(widths)
      rejection_occurred = self._outlier_detector(self._reflections)

    # set the centroid_outlier flag in the original indexed reflections
    ioutliers = self._reflections.get_flags(self._reflections.flags.centroid_outlier)
    ioutliers = self._reflections['iobs'].select(ioutliers)
    self._indexed.set_flags(ioutliers, self._indexed.flags.centroid_outlier)

    msg = "Removing reflections not matched to predictions"
    if rejection_occurred: msg += " or marked as outliers"
    logger.debug(msg)

    # delete all reflections from the manager that do not have a prediction
    # or were flagged as outliers
    has_pred = self._reflections.get_flags(self._reflections.flags.used_in_refinement)
    inlier = ~self._reflections.get_flags(self._reflections.flags.centroid_outlier)
    self._reflections = self._reflections.select(has_pred & inlier)

    logger.debug("%d reflections remain in the manager", len(self._reflections))

    # print summary after outlier rejection
    if rejection_occurred and self._verbosity > 1: self.print_stats_on_matches()

    # form working and free subsets
    self._create_working_set()

    logger.debug("Working set size = %d observations", self.get_sample_size())

    return

  def _id_refs_to_keep(self, obs_data):
    """Create a selection of observations that pass certain conditions.

    This step includes rejection of reflections too close to the spindle,
    reflections measured outside the scan range, rejection of the (0,0,0)
    Miller index and rejection of reflections with the overload flag set.
    Outlier rejection is done later."""

    # first exclude reflections with miller index set to 0,0,0
    sel1 = obs_data['miller_index'] != (0,0,0)

    # exclude reflections with overloads, as these have worse centroids
    sel2 = ~obs_data.get_flags(obs_data.flags.overloaded)

    # combine selections
    sel = sel1 & sel2
    inc = flex.size_t_range(len(obs_data)).select(sel)
    obs_data = obs_data.select(sel)

    # Default to True to pass the following test if there is no rotation axis
    # for a particular experiment
    to_keep = flex.bool(len(inc), True)

    for iexp, exp in enumerate(self._experiments):
      axis = self._axes[iexp]
      if not axis or exp.scan is None: continue
      if exp.scan.get_oscillation()[1] == 0.0: continue
      sel = obs_data['id'] == iexp
      s0 = self._s0vecs[iexp]
      s1 = obs_data['s1'].select(sel)
      phi = obs_data['xyzobs.mm.value'].parts()[2].select(sel)

      # first test: reject reflections for which the parallelepiped formed
      # between the gonio axis, s0 and s1 has a volume of less than the cutoff.
      # Those reflections are by definition closer to the spindle-beam
      # plane and for low values of the cutoff are troublesome to
      # integrate anyway.
      p_vol = flex.abs(s1.cross(flex.vec3_double(s1.size(), s0)).dot(axis))
      passed1 = p_vol > self._close_to_spindle_cutoff

      # second test: reject reflections that lie outside the scan range
      passed2 = exp.scan.is_angle_valid(phi, deg=False)

      # sanity check to catch a mutilated scan that does not make sense
      if passed2.count(True) == 0:
        from libtbx.utils import Sorry
        raise Sorry("Experiment id {0} contains no reflections with valid "
                    "scan angles".format(iexp))

      # combine tests so far
      to_update = passed1 & passed2

      # third test: reject reflections close to the centres of the first and
      # last images in the scan
      if self._trim_scan_edges > 0.0:
        edge1, edge2 = [e + 0.5 for e in exp.scan.get_image_range()]
        edge1 = exp.scan.get_angle_from_image_index(edge1, deg=False)
        edge1 += (self._trim_scan_edges)
        edge2 = exp.scan.get_angle_from_image_index(edge2, deg=False)
        edge2 -= (self._trim_scan_edges)
        passed3 = ((edge1 <= phi) & (phi <= edge2))

        # combine the last test only if there would be a reasonable number of
        # reflections left for refinement
        tmp = to_update
        to_update = to_update & passed3
        if to_update.count(True) < 40:
          logger.warning("Too few reflections to trim centroids from the scan "
            "edges. Resetting trim_scan_edges=0.0")
          to_update = tmp

      # make selection
      to_keep.set_selected(sel, to_update)

    inc = inc.select(to_keep)

    return inc

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

    return len(self._reflections)

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
    nref = len(l)

    from libtbx.table_utils import simple_table
    from scitbx.math import five_number_summary
    x_resid = l['x_resid']
    y_resid = l['y_resid']
    phi_resid = l['phi_resid']
    w_x, w_y, w_phi = l['xyzobs.mm.weights'].parts()

    msg = "\nSummary statistics for {0} observations".format(nref) +\
          " matched to predictions:"
    header = ["", "Min", "Q1", "Med", "Q3", "Max"]
    rows = []
    try:
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
    except IndexError:
      # zero length reflection list
      logger.warning("Unable to calculate summary statistics for zero observations")
      return
    logger.info(msg)
    logger.info(st.format())
    logger.info("")

    # sorting is expensive and the following table is only of interest in
    # special cases, so return now if verbosity is not high
    if self._verbosity < 3: return

    if nref < 20:
      logger.debug("Fewer than 20 reflections matched!")
      return

    sl = self._sort_obs_by_residual(l)
    logger.debug("Reflections with the worst 20 positional residuals:")
    header = ['Miller index', 'x_resid', 'y_resid', 'phi_resid', 'pnl',
              'x_obs', 'y_obs', 'phi_obs', 'x_obs\nweight', 'y_obs\nweight',
              'phi_obs\nweight']
    rows = []
    for i in xrange(20):
      e = sl[i]
      x_obs, y_obs, phi_obs = e['xyzobs.mm.value']
      rows.append(['% 3d, % 3d, % 3d'%e['miller_index'],
                   '%5.3f'%e['x_resid'],
                   '%5.3f'%e['y_resid'],
                   '%6.4f'%(e['phi_resid'] * RAD2DEG),
                   '%d'%e['panel'],
                   '%5.3f'%x_obs,
                   '%5.3f'%y_obs,
                   '%6.4f'%(phi_obs * RAD2DEG),
                   '%5.3f'%e['xyzobs.mm.weights'][0],
                   '%5.3f'%e['xyzobs.mm.weights'][1],
                   '%6.4f'%(e['xyzobs.mm.weights'][2] * DEG2RAD**2)])
    logger.debug(simple_table(rows, header).format())

    sl = self._sort_obs_by_residual(sl, angular=True)
    logger.debug("\nReflections with the worst 20 angular residuals:")
    rows=[]
    for i in xrange(20):
      e = sl[i]
      x_obs, y_obs, phi_obs = e['xyzobs.mm.value']
      rows.append(['% 3d, % 3d, % 3d'%e['miller_index'],
                   '%5.3f'%e['x_resid'],
                   '%5.3f'%e['y_resid'],
                   '%6.4f'%(e['phi_resid'] * RAD2DEG),
                   '%d'%e['panel'],
                   '%5.3f'%x_obs,
                   '%5.3f'%y_obs,
                   '%6.4f'%(phi_obs * RAD2DEG),
                   '%5.3f'%e['xyzobs.mm.weights'][0],
                   '%5.3f'%e['xyzobs.mm.weights'][1],
                   '%6.4f'%(e['xyzobs.mm.weights'][2] * DEG2RAD**2)])
    logger.debug(simple_table(rows, header).format())
    logger.debug("")

    return

  def reset_accepted_reflections(self, reflections=None):
    """Reset use flags for all observations in preparation for a new set of
    predictions"""

    # if not passing in reflections, take the internally managed table
    if reflections is None: reflections = self._reflections

    mask = reflections.get_flags(reflections.flags.used_in_refinement)
    reflections.unset_flags(mask, reflections.flags.used_in_refinement)
    return

  def get_obs(self):
    """Get the list of managed observations"""

    return self._reflections

  def filter_obs(self, sel):
    '''Perform a flex array selection on the managed observations, so that
    external classes can filter according to criteria not available here'''

    self._reflections = self._reflections.select(sel)
    return self._reflections

class StillsReflectionManager(ReflectionManager):
  """Overloads for a Reflection Manager that does not exclude
  reflections too close to the spindle, and reports only information
  about X, Y, DelPsi residuals"""

  _weighting_strategy = weighting_strategies.StillsWeightingStrategy()

  def print_stats_on_matches(self):
    """Print some basic statistics on the matches"""

    l = self.get_matches()
    nref = len(l)

    from libtbx.table_utils import simple_table
    from scitbx.math import five_number_summary
    x_resid = l['x_resid']
    y_resid = l['y_resid']
    delpsi = l['delpsical.rad']
    w_x, w_y, _ = l['xyzobs.mm.weights'].parts()
    w_delpsi = l['delpsical.weights']

    msg = "\nSummary statistics for {0} observations".format(nref) +\
          " matched to predictions:"
    header = ["", "Min", "Q1", "Med", "Q3", "Max"]
    rows = []
    try:
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
    except IndexError:
      # zero length reflection list
      logger.warning("Unable to calculate summary statistics for zero observations")
      return
    logger.info(msg)
    st = simple_table(rows, header)
    logger.info(st.format())
    logger.info("")

    # sorting is expensive and the following table is only of interest in
    # special cases, so return now if verbosity is not high
    if self._verbosity < 3: return

    if nref < 20:
      logger.debug("Fewer than 20 reflections matched!")
      return

    sl = self._sort_obs_by_residual(l)
    logger.debug("Reflections with the worst 20 positional residuals:")
    header = ['Miller index', 'x_resid', 'y_resid', 'pnl',
              'x_obs', 'y_obs', 'x_obs\nweight', 'y_obs\nweight']
    rows = []
    for i in xrange(20):
      e = sl[i]
      x_obs, y_obs, _ = e['xyzobs.mm.value']
      rows.append(['% 3d, % 3d, % 3d'%e['miller_index'],
                   '%5.3f'%e['x_resid'],
                   '%5.3f'%e['y_resid'],
                   '%d'%e['panel'],
                   '%5.3f'%x_obs,
                   '%5.3f'%y_obs,
                   '%5.3f'%e['xyzobs.mm.weights'][0],
                   '%5.3f'%e['xyzobs.mm.weights'][1]])
    logger.debug(simple_table(rows, header).format())
    logger.debug("")

    return
