from __future__ import absolute_import, division, print_function

# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1

import cmath
import copy
import logging
import math
import random
import sys

from libtbx import adopt_init_args
from libtbx import easy_mp
from libtbx.test_utils import approx_equal
from libtbx.utils import plural_s
from scitbx import matrix
from scitbx.simplex import simplex_opt
from scitbx.array_family import flex
import iotbx.phil
from rstbx.indexing_api import dps_extended
from rstbx.indexing_api.lattice import DPS_primitive_lattice
from rstbx.dps_core import Direction, Directional_FFT

from dxtbx.serialize import dump
from dials.algorithms.indexing.indexer import find_max_cell
from dials.util import log
from dials.util import Sorry
from dials.util.options import OptionParser
from dials.util.options import flatten_experiments, flatten_reflections
from dials.util.slice import slice_reflections

logger = logging.getLogger("dials.command_line.search_beam_position")

help_message = """

Search for a better beam centre using the results of spot finding. Based on
method of Sauter et al., J. Appl. Cryst. 37, 399-409 (2004).

Examples::

  dials.search_beam_position imported.expt strong.refl

"""

phil_scope = iotbx.phil.parse(
    """
nproc = Auto
  .type = int(value_min=1)
plot_search_scope = False
  .type = bool
max_cell = None
  .type = float
  .help = "Known max cell (otherwise will compute from spot positions)"
image_range = None
  .help = "The range of images to use in indexing. Number of arguments"
    "must be a factor of two. Specifying \"0 0\" will use all images"
    "by default. The given range follows C conventions"
    "(e.g. j0 <= j < j1)."
  .type = ints(size=2)
  .multiple = True
max_reflections = 10000
  .type = int(value_min=1)
  .help = "Maximum number of reflections to use in the search for better"
          "experimental model. If the number of input reflections is greater"
          "then a random subset of reflections will be used."
mm_search_scope = 4.0
  .help = "Global radius of origin offset search."
  .type = float(value_min=0)
wide_search_binning = 2
  .help = "Modify the coarseness of the wide grid search for the beam centre."
  .type = float(value_min=0)
n_macro_cycles = 1
  .type = int
  .help = "Number of macro cycles for an iterative beam centre search."

seed = 42
  .type = int(value_min=0)

output {
  experiments = optimised.expt
    .type = path
  log = "dials.search_beam_position.log"
    .type = str
  debug_log = "dials.search_beam_position.debug.log"
    .type = str
}
"""
)

master_params = phil_scope.fetch().extract()

dps_phil_scope = iotbx.phil.parse(
    """
include scope rstbx.phil.phil_preferences.indexing_api_defs
d_min = None
  .type = float(value_min=0)
""",
    process_includes=True,
)


class better_experimental_model_discovery(object):
    def __init__(
        self,
        imagesets,
        spot_lists,
        solution_lists,
        amax_lists,
        horizon_phil,
        wide_search_binning=1,
    ):
        adopt_init_args(self, locals())

    def optimize_origin_offset_local_scope(self):
        """Local scope: find the optimal origin-offset closest to the current overall detector position
        (local minimum, simple minimization)"""

        beam = self.imagesets[0].get_beam()
        s0 = matrix.col(beam.get_s0())
        # construct two vectors that are perpendicular to the beam.  Gives a basis for refining beam
        axis = matrix.col((1, 0, 0))
        beamr0 = s0.cross(axis).normalize()
        beamr1 = beamr0.cross(s0).normalize()
        beamr2 = beamr1.cross(s0).normalize()

        assert approx_equal(s0.dot(beamr1), 0.0)
        assert approx_equal(s0.dot(beamr2), 0.0)
        assert approx_equal(beamr2.dot(beamr1), 0.0)
        # so the orthonormal vectors are self.S0_vector, beamr1 and beamr2

        if self.horizon_phil.indexing.mm_search_scope:
            scope = self.horizon_phil.indexing.mm_search_scope
            plot_px_sz = self.imagesets[0].get_detector()[0].get_pixel_size()[0]
            plot_px_sz *= self.wide_search_binning
            grid = max(1, int(scope / plot_px_sz))
            widegrid = 2 * grid + 1
            scores = flex.double()
            for y in range(-grid, grid + 1):
                for x in range(-grid, grid + 1):
                    new_origin_offset = (
                        x * plot_px_sz * beamr1 + y * plot_px_sz * beamr2
                    )
                    score = 0
                    for i in range(len(self.imagesets)):
                        score += self.get_origin_offset_score(
                            new_origin_offset,
                            self.solution_lists[i],
                            self.amax_lists[i],
                            self.spot_lists[i],
                            self.imagesets[i],
                        )
                    scores.append(score)

            def igrid(x):
                return x - (widegrid // 2)

            idxs = [igrid(i) * plot_px_sz for i in range(widegrid)]

            # if there are several similarly high scores, then choose the closest
            # one to the current beam centre
            potential_offsets = flex.vec3_double()
            if scores.all_eq(0):
                raise Sorry("No valid scores")
            sel = scores > (0.9 * flex.max(scores))
            for i in sel.iselection():
                offset = (idxs[i % widegrid]) * beamr1 + (idxs[i // widegrid]) * beamr2
                potential_offsets.append(offset.elems)
                # print offset.length(), scores[i]
            wide_search_offset = matrix.col(
                potential_offsets[flex.min_index(potential_offsets.norms())]
            )

        else:
            wide_search_offset = None

        # DO A SIMPLEX MINIMIZATION

        class test_simplex_method(object):
            def __init__(selfOO, wide_search_offset=None):
                selfOO.starting_simplex = []
                selfOO.n = 2
                selfOO.wide_search_offset = wide_search_offset
                for ii in range(selfOO.n + 1):
                    selfOO.starting_simplex.append(flex.random_double(selfOO.n))
                selfOO.optimizer = simplex_opt(
                    dimension=selfOO.n,
                    matrix=selfOO.starting_simplex,
                    evaluator=selfOO,
                    tolerance=1e-7,
                )
                selfOO.x = selfOO.optimizer.get_solution()
                selfOO.offset = selfOO.x[0] * 0.2 * beamr1 + selfOO.x[1] * 0.2 * beamr2
                if selfOO.wide_search_offset is not None:
                    selfOO.offset += selfOO.wide_search_offset

            def target(selfOO, vector):
                trial_origin_offset = (
                    vector[0] * 0.2 * beamr1 + vector[1] * 0.2 * beamr2
                )
                if selfOO.wide_search_offset is not None:
                    trial_origin_offset += selfOO.wide_search_offset
                target = 0
                for i in range(len(self.imagesets)):
                    target -= self.get_origin_offset_score(
                        trial_origin_offset,
                        self.solution_lists[i],
                        self.amax_lists[i],
                        self.spot_lists[i],
                        self.imagesets[i],
                    )
                return target

        MIN = test_simplex_method(wide_search_offset=wide_search_offset)
        new_offset = MIN.offset

        if self.horizon_phil.indexing.plot_search_scope:
            scope = self.horizon_phil.indexing.mm_search_scope
            plot_px_sz = self.imagesets[0].get_detector()[0].get_pixel_size()[0]
            grid = max(1, int(scope / plot_px_sz))
            scores = flex.double()
            for y in range(-grid, grid + 1):
                for x in range(-grid, grid + 1):
                    new_origin_offset = (
                        x * plot_px_sz * beamr1 + y * plot_px_sz * beamr2
                    )
                    score = 0
                    for i in range(len(self.imagesets)):
                        score += self.get_origin_offset_score(
                            new_origin_offset,
                            self.solution_lists[i],
                            self.amax_lists[i],
                            self.spot_lists[i],
                            self.imagesets[i],
                        )
                    scores.append(score)

            def show_plot(widegrid, excursi):
                excursi.reshape(flex.grid(widegrid, widegrid))
                idx_max = flex.max_index(excursi)

                def igrid(x):
                    return x - (widegrid // 2)

                idxs = [igrid(i) * plot_px_sz for i in range(widegrid)]

                from matplotlib import pyplot as plt

                plt.figure()
                CS = plt.contour(
                    [igrid(i) * plot_px_sz for i in range(widegrid)],
                    [igrid(i) * plot_px_sz for i in range(widegrid)],
                    excursi.as_numpy_array(),
                )
                plt.clabel(CS, inline=1, fontsize=10, fmt="%6.3f")
                plt.title("Wide scope search for detector origin offset")
                plt.scatter([0.0], [0.0], color="g", marker="o")
                plt.scatter([new_offset[0]], [new_offset[1]], color="r", marker="*")
                plt.scatter(
                    [idxs[idx_max % widegrid]],
                    [idxs[idx_max // widegrid]],
                    color="k",
                    marker="s",
                )
                plt.axes().set_aspect("equal")
                plt.xlabel("offset (mm) along beamr1 vector")
                plt.ylabel("offset (mm) along beamr2 vector")
                plt.savefig("search_scope.png")

                # changing value
                trial_origin_offset = (idxs[idx_max % widegrid]) * beamr1 + (
                    idxs[idx_max // widegrid]
                ) * beamr2
                return trial_origin_offset

            show_plot(widegrid=2 * grid + 1, excursi=scores)

        return dps_extended.get_new_detector(
            self.imagesets[0].get_detector(), new_offset
        )

    def get_origin_offset_score(
        self, trial_origin_offset, solutions, amax, spots_mm, imageset
    ):
        trial_detector = dps_extended.get_new_detector(
            imageset.get_detector(), trial_origin_offset
        )

        # Key point for this is that the spots must correspond to detector
        # positions not to the correct RS position => reset any fixed rotation
        # to identity - copy in case called from elsewhere

        gonio = copy.deepcopy(imageset.get_goniometer())
        gonio.set_fixed_rotation((1, 0, 0, 0, 1, 0, 0, 0, 1))
        spots_mm.map_centroids_to_reciprocal_space(
            trial_detector, imageset.get_beam(), gonio
        )

        return self.sum_score_detail(spots_mm["rlp"], solutions, amax=amax)

    def sum_score_detail(
        self, reciprocal_space_vectors, solutions, granularity=None, amax=None
    ):
        """Evaluates the probability that the trial value of (S0_vector | origin_offset) is correct,
        given the current estimate and the observations.  The trial value comes through the
        reciprocal space vectors, and the current estimate comes through the short list of
        DPS solutions. Actual return value is a sum of NH terms, one for each DPS solution, each ranging
        from -1.0 to 1.0"""

        nh = min(solutions.size(), 20)  # extended API
        sum_score = 0.0
        for t in range(nh):
            # if t!=unique:continue
            dfft = Directional_FFT(
                angle=Direction(solutions[t]),
                xyzdata=reciprocal_space_vectors,
                granularity=5.0,
                amax=amax,  # extended API XXX These values have to come from somewhere!
                F0_cutoff=11,
            )
            kval = dfft.kval()
            kmax = dfft.kmax()
            # kval_cutoff = self.raw_spot_input.size()/4.0; # deprecate record
            kval_cutoff = reciprocal_space_vectors.size() / 4.0
            # deprecate record
            if kval > kval_cutoff:
                ff = dfft.fft_result
                kbeam = ((-dfft.pmin) / dfft.delta_p) + 0.5
                Tkmax = cmath.phase(ff[kmax])
                backmax = math.cos(
                    Tkmax + (2 * math.pi * kmax * kbeam / (2 * ff.size() - 1))
                )
                ### Here it should be possible to calculate a gradient.
                ### Then minimize with respect to two coordinates.  Use lbfgs?  Have second derivatives?
                ### can I do something local to model the cosine wave?
                ### direction of wave travel.  Period. phase.
                sum_score += backmax
            # if t == unique:
            #  print t, kmax, dfft.pmin, dfft.delta_p, Tkmax,(2*math.pi*kmax*kbeam/(2*ff.size()-1))
        return sum_score


def run_dps(args):
    imageset, spots_mm, max_cell, params = args

    detector = imageset.get_detector()
    beam = imageset.get_beam()
    goniometer = imageset.get_goniometer()

    # max_cell: max possible cell in Angstroms; set to None, determine from data
    # recommended_grid_sampling_rad: grid sampling in radians; guess for now

    DPS = DPS_primitive_lattice(
        max_cell=max_cell, recommended_grid_sampling_rad=None, horizon_phil=params
    )

    DPS.S0_vector = matrix.col(beam.get_s0())
    DPS.inv_wave = 1.0 / beam.get_wavelength()
    if goniometer is None:
        DPS.axis = matrix.col((1, 0, 0))
    else:
        DPS.axis = matrix.col(goniometer.get_rotation_axis())
    DPS.set_detector(detector)

    # transform input into what Nick needs
    # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees

    data = flex.vec3_double()
    for spot in spots_mm:
        data.append(
            (
                spot["xyzobs.mm.value"][0],
                spot["xyzobs.mm.value"][1],
                spot["xyzobs.mm.value"][2] * 180.0 / math.pi,
            )
        )

    logger.info("Running DPS using %i reflections" % len(data))

    DPS.index(
        raw_spot_input=data, panel_addresses=flex.int([s["panel"] for s in spots_mm])
    )
    solutions = DPS.getSolutions()

    logger.info(
        "Found %i solution%s with max unit cell %.2f Angstroms."
        % (len(solutions), plural_s(len(solutions))[1], DPS.amax)
    )
    if len(solutions) < 3:

        raise Sorry(
            "Not enough solutions: found %i, need at least 3" % (len(solutions))
        )
    return dict(solutions=flex.vec3_double([s.dvec for s in solutions]), amax=DPS.amax)


def discover_better_experimental_model(
    imagesets, spot_lists, params, dps_params, nproc=1, wide_search_binning=1
):
    assert len(imagesets) == len(spot_lists)
    assert len(imagesets) > 0
    # XXX should check that all the detector and beam objects are the same

    spot_lists_mm = []
    max_cell_list = []

    detector = imagesets[0].get_detector()
    beam = imagesets[0].get_beam()

    beam_panel = detector.get_panel_intersection(beam.get_s0())

    if beam_panel == -1:

        raise Sorry("input beam does not intersect detector")

    for imageset, spots in zip(imagesets, spot_lists):
        if "imageset_id" not in spots:
            spots["imageset_id"] = spots["id"]

        spots_mm = copy.deepcopy(spots)
        spots_mm.centroid_px_to_mm(imageset.get_detector(), scan=imageset.get_scan())

        spots_mm.map_centroids_to_reciprocal_space(
            detector=imageset.get_detector(),
            beam=imageset.get_beam(),
            goniometer=imageset.get_goniometer(),
        )

        if dps_params.d_min is not None:
            d_spacings = 1 / spots_mm["rlp"].norms()
            sel = d_spacings > dps_params.d_min
            spots_mm = spots_mm.select(sel)

        # derive a max_cell from mm spots

        if params.max_cell is None:
            max_cell = find_max_cell(
                spots_mm, max_cell_multiplier=1.3, step_size=45
            ).max_cell
            max_cell_list.append(max_cell)

        if (
            params.max_reflections is not None
            and spots_mm.size() > params.max_reflections
        ):
            logger.info(
                "Selecting subset of %i reflections for analysis"
                % params.max_reflections
            )
            perm = flex.random_permutation(spots_mm.size())
            sel = perm[: params.max_reflections]
            spots_mm = spots_mm.select(sel)

        spot_lists_mm.append(spots_mm)

    if params.max_cell is None:
        max_cell = flex.median(flex.double(max_cell_list))
    else:
        max_cell = params.max_cell
    args = [
        (imageset, spots, max_cell, dps_params)
        for imageset, spots in zip(imagesets, spot_lists_mm)
    ]

    results = easy_mp.parallel_map(
        func=run_dps,
        iterable=args,
        processes=nproc,
        method="multiprocessing",
        preserve_order=True,
        asynchronous=True,
        preserve_exception_message=True,
    )
    solution_lists = [r["solutions"] for r in results]
    amax_list = [r["amax"] for r in results]
    assert len(solution_lists) > 0

    detector = imagesets[0].get_detector()
    beam = imagesets[0].get_beam()

    # perform calculation
    if dps_params.indexing.improve_local_scope == "origin_offset":
        discoverer = better_experimental_model_discovery(
            imagesets,
            spot_lists_mm,
            solution_lists,
            amax_list,
            dps_params,
            wide_search_binning=wide_search_binning,
        )
        new_detector = discoverer.optimize_origin_offset_local_scope()
        old_panel, old_beam_centre = detector.get_ray_intersection(beam.get_s0())
        new_panel, new_beam_centre = new_detector.get_ray_intersection(beam.get_s0())

        old_beam_centre_px = detector[old_panel].millimeter_to_pixel(old_beam_centre)
        new_beam_centre_px = new_detector[new_panel].millimeter_to_pixel(
            new_beam_centre
        )

        logger.info(
            "Old beam centre: %.2f, %.2f mm" % old_beam_centre
            + " (%.1f, %.1f px)" % old_beam_centre_px
        )
        logger.info(
            "New beam centre: %.2f, %.2f mm" % new_beam_centre
            + " (%.1f, %.1f px)" % new_beam_centre_px
        )
        logger.info(
            "Shift: %.2f, %.2f mm"
            % (matrix.col(old_beam_centre) - matrix.col(new_beam_centre)).elems
            + " (%.1f, %.1f px)"
            % (matrix.col(old_beam_centre_px) - matrix.col(new_beam_centre_px)).elems
        )
        return new_detector, beam
    elif dps_params.indexing.improve_local_scope == "S0_vector":
        raise NotImplementedError()


def run(args):
    usage = "dials.search_beam_position [options] imported.expt strong.refl"

    parser = OptionParser(
        usage=usage,
        phil=phil_scope,
        read_experiments=True,
        read_reflections=True,
        check_format=False,
        epilog=help_message,
    )

    params, options = parser.parse_args(show_diff_phil=False)
    experiments = flatten_experiments(params.input.experiments)
    reflections = flatten_reflections(params.input.reflections)

    if len(experiments) == 0 or len(reflections) == 0:
        parser.print_help()
        exit(0)

    # Configure the logging
    log.config(info=params.output.log, debug=params.output.debug_log)

    # Log the diff phil
    diff_phil = parser.diff_phil.as_str()
    if diff_phil != "":
        logger.info("The following parameters have been modified:\n")
        logger.info(diff_phil)

    if params.seed is not None:
        flex.set_random_seed(params.seed)
        random.seed(params.seed)

    imagesets = experiments.imagesets()
    # Split all the refln tables by ID, corresponding to the respective imagesets
    reflections = [
        refl_unique_id
        for refl in reflections
        for refl_unique_id in refl.split_by_experiment_id()
    ]

    assert len(imagesets) > 0
    assert len(reflections) == len(imagesets)

    if params.image_range is not None and len(params.image_range) > 0:
        reflections = [
            slice_reflections(refl, params.image_range) for refl in reflections
        ]

    dps_params = dps_phil_scope.extract()
    # for development, we want an exhaustive plot of beam probability map:
    dps_params.indexing.plot_search_scope = params.plot_search_scope
    dps_params.indexing.mm_search_scope = params.mm_search_scope

    for i in range(params.n_macro_cycles):
        if params.n_macro_cycles > 1:
            logger.info("Starting macro cycle %i" % (i + 1))
        new_detector, new_beam = discover_better_experimental_model(
            imagesets,
            reflections,
            params,
            dps_params,
            nproc=params.nproc,
            wide_search_binning=params.wide_search_binning,
        )
        for experiment in experiments:
            experiment.beam = new_beam
            experiment.detector = new_detector
        logger.info("")

    logger.info("Saving optimized experiments to %s" % params.output.experiments)
    dump.experiment_list(experiments, params.output.experiments)


if __name__ == "__main__":
    run(sys.argv[1:])
