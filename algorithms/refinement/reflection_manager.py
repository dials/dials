"""Contains classes used to manage the reflections used during refinement,
principally ReflectionManager."""
from __future__ import absolute_import, division, print_function

import copy
import logging
import math
import random

from dials.algorithms.refinement import DialsRefineConfigError
from dials.algorithms.refinement import weighting_strategies
from dials.algorithms.refinement.analysis.centroid_analysis import CentroidAnalyser
from dials.algorithms.refinement.outlier_detection.outlier_base import (
    phil_str as outlier_phil_str,
)
from dials.algorithms.refinement.refinement_helpers import (
    calculate_frame_numbers,
    set_obs_s1,
)
from dials.array_family import flex
import libtbx
from libtbx.phil import parse
from libtbx.table_utils import simple_table
from scitbx import matrix
from scitbx.math import five_number_summary

logger = logging.getLogger(__name__)

# constants
RAD2DEG = 180.0 / math.pi
DEG2RAD = math.pi / 180.0

# PHIL
format_data = {"outlier_phil": outlier_phil_str}
phil_str = (
    """
    reflections_per_degree = Auto
      .help = "The number of centroids per degree of the sweep to use in"
              "refinement. The default (Auto) uses all reflections unless"
              "the dataset is wider than a single turn. Then the number of"
              "reflections may be reduced until a minimum of 100 per degree of"
              "the sweep is reached to speed up calculations. Set this to None"
              "to force use all of suitable reflections."
      .type = float(value_min=0.)
      .expert_level = 1

    minimum_sample_size = 1000
      .help = "cutoff that determines whether subsetting of the input"
              "reflection list is done"
      .type = int
      .expert_level = 1

    maximum_sample_size = None
      .help = "The maximum number of reflections to use in refinement."
              "Overrides reflections_per_degree if that produces a"
              "larger sample size."
      .type = int(value_min=1)
      .expert_level = 1

    random_seed = 42
      .help = "Random seed to use when sampling to create a working set of"
              "reflections. May be int or None."
      .type = int
      .expert_level = 1

    close_to_spindle_cutoff = 0.02
      .help = "The inclusion criterion currently uses the volume of the"
              "parallelepiped formed by the spindle axis, the incident"
              "beam and the scattered beam. If this is lower than some"
              "value then the reflection is excluded from refinement."
              "In detector space, these are the reflections located close"
              "to the rotation axis."
      .type = float(value_min = 0)
      .expert_level = 1

    trim_scan_edges = 0.0
      .help = "Reflections within this value in degrees from the centre of the"
              "first or last image of the scan will be removed before"
              "refinement, unless doing so would result in too few remaining"
              "reflections. Reflections that are truncated at the scan edges"
              "have poorly-determined centroids and can bias the refined model"
              "if they are included."
      .type = float(value_min=0,value_max=1)
      .expert_level = 1

    weighting_strategy
      .help = "Parameters to configure weighting strategy overrides"
      .expert_level = 1
    {
      override = statistical stills constant external_deltapsi
        .help = "selection of a strategy to override default weighting behaviour"
        .type = choice

      delpsi_constant = 1000000
        .help = "used by the stills strategy to choose absolute weight value"
                "for the angular distance from Ewald sphere term of the target"
                "function, whilst the X and Y parts use statistical weights"
        .type = float(value_min = 0)

      constants = 1.0 1.0 1.0
        .help = "constant weights for three parts of the target function,"
                "whether the case is for stills or scans. The default gives"
                "unit weighting."
        .type = floats(size = 3, value_min = 0)
    }

    %(outlier_phil)s
"""
    % format_data
)
phil_scope = parse(phil_str)


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

    def _create_block_columns(self):
        """Create a column to contain the block number."""

        from scitbx.array_family import flex

        self._reflections["block"] = flex.size_t(len(self._reflections))
        self._reflections["block_centre"] = flex.double(len(self._reflections))

    @staticmethod
    def _check_scan_range(exp_phi, scan):
        """Check that the observed reflections fill the scan-range"""

        # Allow up to 5 degrees between the observed phi extrema and the
        # scan edges
        start, stop = scan.get_oscillation_range(deg=False)
        if min(exp_phi) - start > 0.087266 or stop - max(exp_phi) > 0.087266:
            raise DialsRefineConfigError("The reflections do not fill the scan range.")

    def per_width(self, width, deg=True):
        """Set blocks for all experiments according to a constant width"""

        if deg:
            width *= DEG2RAD
        self._create_block_columns()

        # get observed phi in radians
        phi_obs = self._reflections["xyzobs.mm.value"].parts()[2]

        for iexp, exp in enumerate(self._experiments):

            sel = self._reflections["id"] == iexp
            isel = sel.iselection()
            exp_phi = phi_obs.select(isel)
            self._check_scan_range(exp_phi, exp.scan)

            start, stop = exp.scan.get_oscillation_range(deg=False)
            nblocks = int(abs(stop - start) / width) + 1
            # ensure width has the right sign and is wide enough that all reflections
            # get assigned a block
            _width = cmp(stop, start) * width + 1e-11
            half_width = width * (0.5 - 1e-11)  # ensure round down behaviour

            block_starts = [start + n * _width for n in range(nblocks)]
            block_centres = [
                exp.scan.get_array_index_from_angle(e + half_width, deg=False)
                for e in block_starts
            ]

            for b_num, (b_start, b_cent) in enumerate(zip(block_starts, block_centres)):
                sub_isel = isel.select(
                    (b_start <= exp_phi) & (exp_phi <= (b_start + _width))
                )
                self._reflections["block"].set_selected(sub_isel, b_num)
                self._reflections["block_centre"].set_selected(sub_isel, b_cent)

        return self._reflections

    def per_image(self):
        """Set one block per image for all experiments"""

        self._create_block_columns()

        # get observed phi in radians
        phi_obs = self._reflections["xyzobs.mm.value"].parts()[2]

        for iexp, exp in enumerate(self._experiments):

            sel = self._reflections["id"] == iexp
            isel = sel.iselection()
            exp_phi = phi_obs.select(isel)
            self._check_scan_range(exp_phi, exp.scan)

            # convert phi to integer frames
            frames = exp.scan.get_array_index_from_angle(exp_phi, deg=False)
            frames = flex.floor(frames).iround()

            start, stop = flex.min(frames), flex.max(frames)
            frame_range = range(start, stop + 1)

            for f_num, f in enumerate(frame_range):
                sub_isel = isel.select(frames == f)
                f_cent = f + 0.5
                self._reflections["block"].set_selected(sub_isel, f_num)
                self._reflections["block_centre"].set_selected(sub_isel, f_cent)

        return self._reflections


class ReflectionManagerFactory(object):
    @staticmethod
    def from_parameters_reflections_experiments(
        params, reflections, experiments, do_stills=False
    ):

        """Given a set of parameters and models, build a reflection manager

        Params:
            params The input parameters

        Returns:
            The reflection manager instance
        """

        # While a random subset of reflections is used, continue to
        # set random.seed to get consistent behaviour
        if params.random_seed is not None:
            random.seed(params.random_seed)
            flex.set_random_seed(params.random_seed)
            logger.debug("Random seed set to %d", params.random_seed)

        # check whether we deal with stills or scans
        if do_stills:
            refman = StillsReflectionManager
            # check incompatible weighting strategy
            if params.weighting_strategy.override == "statistical":
                raise DialsRefineConfigError(
                    'The "statistical" weighting strategy is not compatible '
                    "with stills refinement"
                )
        else:
            refman = ReflectionManager
            # check incompatible weighting strategy
            if params.weighting_strategy.override in ["stills", "external_deltapsi"]:
                msg = (
                    'The "{0}" weighting strategy is not compatible with '
                    "scan refinement"
                ).format(params.weighting_strategy.override)
                raise DialsRefineConfigError(msg)

        # set automatic outlier rejection options
        if params.outlier.algorithm in ("auto", libtbx.Auto):
            if do_stills:
                params.outlier.algorithm = "sauter_poon"
            else:
                params.outlier.algorithm = "mcd"

        if params.outlier.separate_panels is libtbx.Auto:
            if do_stills:
                params.outlier.separate_panels = False
            else:
                params.outlier.separate_panels = True

        if params.outlier.algorithm == "sauter_poon":
            if params.outlier.sauter_poon.px_sz is libtbx.Auto:
                # get this from the first panel of the first detector
                params.outlier.sauter_poon.px_sz = experiments.detectors()[0][
                    0
                ].get_pixel_size()

        # do outlier rejection?
        if params.outlier.algorithm in ("null", None):
            outlier_detector = None
        else:
            if do_stills:
                colnames = ["x_resid", "y_resid"]
                params.outlier.block_width = None
            else:
                colnames = ["x_resid", "y_resid", "phi_resid"]
            from dials.algorithms.refinement.outlier_detection import (
                CentroidOutlierFactory,
            )

            outlier_detector = CentroidOutlierFactory.from_parameters_and_colnames(
                params, colnames
            )

        # override default weighting strategy?
        weighting_strategy = None
        if params.weighting_strategy.override == "statistical":
            from dials.algorithms.refinement.weighting_strategies import (
                StatisticalWeightingStrategy,
            )

            weighting_strategy = StatisticalWeightingStrategy()
        elif params.weighting_strategy.override == "stills":
            from dials.algorithms.refinement.weighting_strategies import (
                StillsWeightingStrategy,
            )

            weighting_strategy = StillsWeightingStrategy(
                params.weighting_strategy.delpsi_constant
            )
        elif params.weighting_strategy.override == "external_deltapsi":
            from dials.algorithms.refinement.weighting_strategies import (
                ExternalDelPsiWeightingStrategy,
            )

            weighting_strategy = ExternalDelPsiWeightingStrategy()
        elif params.weighting_strategy.override == "constant":
            from dials.algorithms.refinement.weighting_strategies import (
                ConstantWeightingStrategy,
            )

            weighting_strategy = ConstantWeightingStrategy(
                *params.weighting_strategy.constants, stills=do_stills
            )

        return refman(
            reflections=reflections,
            experiments=experiments,
            nref_per_degree=params.reflections_per_degree,
            max_sample_size=params.maximum_sample_size,
            min_sample_size=params.minimum_sample_size,
            close_to_spindle_cutoff=params.close_to_spindle_cutoff,
            trim_scan_edges=params.trim_scan_edges,
            outlier_detector=outlier_detector,
            weighting_strategy_override=weighting_strategy,
        )


class ReflectionManager(object):
    """A class to maintain information about observed and predicted
    reflections for refinement.

    This new version keeps the reflections as a reflection table. Initialisation
    is not complete until the ReflectionManager is paired with a target function.
    Then, prediction can be done, followed by outlier rejection and any random
    sampling to form the working subset."""

    _weighting_strategy = weighting_strategies.StatisticalWeightingStrategy()
    experiment_type = "scans"

    def __init__(
        self,
        reflections,
        experiments,
        nref_per_degree=None,
        max_sample_size=None,
        min_sample_size=0,
        close_to_spindle_cutoff=0.02,
        trim_scan_edges=0.0,
        outlier_detector=None,
        weighting_strategy_override=None,
    ):

        if len(reflections) == 0:
            raise ValueError("Empty reflections table provided to ReflectionManager")

        # keep track of models
        self._experiments = experiments
        goniometers = [e.goniometer for e in self._experiments]
        self._axes = [
            matrix.col(g.get_rotation_axis()) if g else None for g in goniometers
        ]
        self._s0vecs = [matrix.col(e.beam.get_s0()) for e in self._experiments]

        # unset the refinement flags (creates flags field if needed)
        reflections.unset_flags(
            flex.size_t_range(len(reflections)),
            flex.reflection_table.flags.used_in_refinement,
        )

        # check that the observed beam vectors are stored: if not, compute them
        n_s1_set = set_obs_s1(reflections, experiments)
        if n_s1_set > 0:
            logger.debug("Set scattering vectors for %d reflections", n_s1_set)

        # keep track of the original indices of the reflections
        reflections["iobs"] = flex.size_t_range(len(reflections))

        # Check for monotonically increasing value range. If not, ref_table isn't sorted,
        # and proceed to sort by id and panel. This is required for the C++ extension
        # modules to allow for nlogn subselection of values used in refinement.
        l_id = reflections["id"]
        id0 = l_id[0]
        for id_x in l_id[1:]:
            if id0 <= id_x:
                id0 = id_x
            else:
                reflections.sort("id")  # Ensuring the ref_table is sorted by id
                reflections.subsort(
                    "id", "panel"
                )  # Ensuring that within each sorted id block, sorting is next performed by panel
                break

        # set up the reflection inclusion criteria
        self._close_to_spindle_cutoff = close_to_spindle_cutoff  # close to spindle
        self._trim_scan_edges = DEG2RAD * trim_scan_edges  # close to the scan edge
        self._outlier_detector = outlier_detector  # for outlier rejection
        self._nref_per_degree = nref_per_degree  # random subsets
        self._max_sample_size = max_sample_size  # sample size ceiling
        self._min_sample_size = min_sample_size  # sample size floor

        # exclude reflections that fail some inclusion criteria
        refs_to_keep = self._id_refs_to_keep(reflections)
        self._accepted_refs_size = len(refs_to_keep)

        # set entering flags for all reflections
        reflections.calculate_entering_flags(self._experiments)

        # set observed frame numbers for all reflections if not already present
        calculate_frame_numbers(reflections, self._experiments)

        # reset all use flags
        self.reset_accepted_reflections(reflections)

        # put full list of indexed reflections aside and select only the reflections
        # that were not excluded to manage
        self._indexed = reflections
        self._reflections = reflections.select(refs_to_keep)

        # set exclusion flag for reflections that failed the tests
        refs_to_excl = flex.bool(len(self._indexed), True)
        refs_to_excl.set_selected(refs_to_keep, False)
        self._indexed.set_flags(
            refs_to_excl, self._indexed.flags.excluded_for_refinement
        )

        # set weights for all kept reflections
        if weighting_strategy_override is not None:
            self._weighting_strategy = weighting_strategy_override
        self._weighting_strategy.calculate_weights(self._reflections)

        # not known until the manager is finalised
        self._sample_size = None

    def get_centroid_analyser(self, debug=False):
        """Create a CentroidAnalysis object for the current reflections"""

        return CentroidAnalyser(self._reflections, debug=debug)

    def finalise(self, analysis=None):
        """Complete initialisation by performing outlier rejection and any
        requested subsetting. If a list of results from a CentroidAnalysis
        object is provided, these may be used to determine outlier rejection
        block widths"""

        logger.debug("Finalising the Reflection Manager")

        # Initially, assume all reflections with predictions can be used
        mask = self._reflections.get_flags(self._reflections.flags.predicted)
        self._reflections.set_flags(mask, self._reflections.flags.used_in_refinement)

        # print summary before outlier rejection
        self.print_stats_on_matches()

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
                    widths = [e.get("block_size") for e in analysis]
                    widths = [max(e, 18.0) if e is not None else None for e in widths]
                    self._outlier_detector.set_block_width(widths)
            rejection_occurred = self._outlier_detector(self._reflections)

        # set the centroid_outlier flag in the original indexed reflections
        ioutliers = self._reflections.get_flags(
            self._reflections.flags.centroid_outlier
        )
        ioutliers = self._reflections["iobs"].select(ioutliers)
        self._indexed.sort("iobs")  # re-sort the indexed reflections
        self._indexed.set_flags(ioutliers, self._indexed.flags.centroid_outlier)

        msg = "Removing reflections not matched to predictions"
        if rejection_occurred:
            msg += " or marked as outliers"
        logger.debug(msg)

        # delete all reflections from the manager that do not have a prediction
        # or were flagged as outliers
        has_pred = self._reflections.get_flags(self._reflections.flags.predicted)
        inlier = ~self._reflections.get_flags(self._reflections.flags.centroid_outlier)
        self._reflections = self._reflections.select(has_pred & inlier)
        self._reflections.set_flags(
            flex.bool(len(self._reflections), True),
            self._reflections.flags.used_in_refinement,
        )

        logger.debug("%d reflections remain in the manager", len(self._reflections))

        # print summary after outlier rejection
        if rejection_occurred:
            self.print_stats_on_matches()

        # form working and free subsets
        self._create_working_set()

        logger.debug("Working set size = %d observations", self.get_sample_size())

    def _id_refs_to_keep(self, obs_data):
        """Create a selection of observations that pass certain conditions.

        This step includes rejection of reflections too close to the spindle,
        reflections measured outside the scan range, rejection of the (0,0,0)
        Miller index and rejection of reflections with the overload flag set.
        Outlier rejection is done later."""

        # first exclude reflections with miller index set to 0,0,0
        sel1 = obs_data["miller_index"] != (0, 0, 0)

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
            if not axis or exp.scan is None:
                continue
            if exp.scan.get_oscillation()[1] == 0.0:
                continue
            sel = obs_data["id"] == iexp
            s0 = self._s0vecs[iexp]
            s1 = obs_data["s1"].select(sel)
            phi = obs_data["xyzobs.mm.value"].parts()[2].select(sel)

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
                raise DialsRefineConfigError(
                    "Experiment id {} contains no reflections with valid "
                    "scan angles".format(iexp)
                )

            # combine tests so far
            to_update = passed1 & passed2

            # third test: reject reflections close to the centres of the first and
            # last images in the scan
            if self._trim_scan_edges > 0.0:
                edge1, edge2 = [e + 0.5 for e in exp.scan.get_image_range()]
                edge1 = exp.scan.get_angle_from_image_index(edge1, deg=False)
                edge1 += self._trim_scan_edges
                edge2 = exp.scan.get_angle_from_image_index(edge2, deg=False)
                edge2 -= self._trim_scan_edges
                passed3 = (edge1 <= phi) & (phi <= edge2)

                # combine the last test only if there would be a reasonable number of
                # reflections left for refinement
                tmp = to_update
                to_update = to_update & passed3
                if to_update.count(True) < 40:
                    logger.warning(
                        "Too few reflections to trim centroids from the scan "
                        "edges. Resetting trim_scan_edges=0.0"
                    )
                    to_update = tmp

            # make selection
            to_keep.set_selected(sel, to_update)

        inc = inc.select(to_keep)

        return inc

    def _create_working_set(self):
        """Make a subset of the indices of reflections to use in refinement"""

        working_isel = flex.size_t()
        for iexp, exp in enumerate(self._experiments):

            sel = self._reflections["id"] == iexp
            isel = sel.iselection()
            # refs = self._reflections.select(sel)
            nrefs = sample_size = len(isel)

            # set sample size according to nref_per_degree (per experiment)
            if exp.scan and self._nref_per_degree:
                sweep_range_rad = exp.scan.get_oscillation_range(deg=False)
                width = abs(sweep_range_rad[1] - sweep_range_rad[0]) * RAD2DEG
                if self._nref_per_degree is libtbx.Auto:
                    # For multi-turn, set sample size to the greater of the approx nref
                    # in a single turn and 100 reflections per degree
                    turns = width / 360.0
                    if turns > 1:
                        approx_nref_1_turn = int(math.ceil(nrefs / turns))
                        sample_size = int(max(approx_nref_1_turn, 100.0 * width))
                else:
                    sample_size = int(self._nref_per_degree * width)

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
        sort_obs = copy.deepcopy(obs)
        if angular:
            sort_obs.sort("phi_resid", reverse=True)
        else:
            sort_obs["key"] = sort_obs["x_resid"] ** 2 + sort_obs["y_resid"] ** 2
            sort_obs.sort("key", reverse=True)
            del sort_obs["key"]
        return sort_obs

    def get_indexed(self):
        """Return the reflections passed in as input"""

        return self._indexed

    def get_matches(self):
        """For every observation used in refinement return (a copy of) all data"""

        return self._reflections.select(
            self._reflections.get_flags(self._reflections.flags.used_in_refinement)
        )

    def get_free_reflections(self):
        """Return all reflections that were accepted for refinement but not chosen
        in the working set"""

        return self._free_reflections

    def print_stats_on_matches(self):
        """Print some basic statistics on the matches"""

        l = self.get_matches()
        nref = len(l)
        if nref == 0:
            logger.warning(
                "Unable to calculate summary statistics for zero observations"
            )
            return

        try:
            x_resid = l["x_resid"]
            y_resid = l["y_resid"]
            phi_resid = l["phi_resid"]
            w_x, w_y, w_phi = l["xyzobs.mm.weights"].parts()
        except RuntimeError:  # wish it were KeyError
            return

        msg = (
            "\nSummary statistics for {} observations".format(nref)
            + " matched to predictions:"
        )
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
        rows.append(["Phi weights"] + ["%.4g" % (e * DEG2RAD ** 2) for e in row_data])
        st = simple_table(rows, header)

        logger.info(msg)
        logger.info(st.format())
        logger.info("")

    def reset_accepted_reflections(self, reflections=None):
        """Reset use flags for all observations in preparation for a new set of
        predictions"""

        # if not passing in reflections, take the internally managed table
        if reflections is None:
            reflections = self._reflections

        mask = reflections.get_flags(reflections.flags.used_in_refinement)
        reflections.unset_flags(mask, reflections.flags.used_in_refinement)

    def get_obs(self):
        """Get the list of managed observations"""

        return self._reflections

    def filter_obs(self, sel):
        """Perform a flex array selection on the managed observations, so that
        external classes can filter according to criteria not available here"""

        self._reflections = self._reflections.select(sel)
        return self._reflections


class StillsReflectionManager(ReflectionManager):
    """Overloads for a Reflection Manager that does not exclude
    reflections too close to the spindle, and reports only information
    about X, Y, DelPsi residuals"""

    _weighting_strategy = weighting_strategies.StillsWeightingStrategy()
    experiment_type = "stills"

    def _id_refs_to_keep(self, obs_data):
        """Create a selection of observations that pass certain conditions.

        Stills-specific version removes checks relevant only to experiments
        with a rotation axis."""

        # first exclude reflections with miller index set to 0,0,0
        sel1 = obs_data["miller_index"] != (0, 0, 0)

        # exclude reflections with overloads, as these have worse centroids
        sel2 = ~obs_data.get_flags(obs_data.flags.overloaded)

        # combine selections
        sel = sel1 & sel2
        inc = flex.size_t_range(len(obs_data)).select(sel)

        return inc

    def print_stats_on_matches(self):
        """Print some basic statistics on the matches"""

        l = self.get_matches()
        nref = len(l)
        if nref == 0:
            logger.warning(
                "Unable to calculate summary statistics for zero observations"
            )
            return

        from libtbx.table_utils import simple_table
        from scitbx.math import five_number_summary

        try:
            x_resid = l["x_resid"]
            y_resid = l["y_resid"]
            delpsi = l["delpsical.rad"]
            w_x, w_y, _ = l["xyzobs.mm.weights"].parts()
            w_delpsi = l["delpsical.weights"]
        except RuntimeError:  # wish it were KeyError
            return

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
        rows.append(
            ["DeltaPsi weights"] + ["%.4g" % (e * DEG2RAD ** 2) for e in row_data]
        )

        msg = (
            "\nSummary statistics for {} observations".format(nref)
            + " matched to predictions:"
        )
        logger.info(msg)
        st = simple_table(rows, header)
        logger.info(st.format())
        logger.info("")
