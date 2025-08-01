from __future__ import annotations

import cmath
import concurrent.futures
import copy
import itertools
import logging
import math

import iotbx.phil
from libtbx.test_utils import approx_equal
from libtbx.utils import plural_s
from rstbx.dps_core import Direction, Directional_FFT
from rstbx.indexing_api import dps_extended
from rstbx.indexing_api.lattice import DPS_primitive_lattice
from rstbx.phil.phil_preferences import indexing_api_defs
from scitbx import matrix
from scitbx.array_family import flex
from scitbx.simplex import simplex_opt

from dials.algorithms.indexing.indexer import find_max_cell
from dials.util import Sorry

logger = logging.getLogger("dials.command_line.search_beam_position")


def optimize_origin_offset_local_scope(
    experiments,
    reflection_lists,
    solution_lists,
    amax_lists,
    mm_search_scope=4,
    wide_search_binning=1,
    plot_search_scope=False,
):
    """Local scope: find the optimal origin-offset closest to the current
    overall detector position (local minimum, simple minimization)"""

    beam = experiments[0].beam
    s0 = matrix.col(beam.get_s0())
    # construct two vectors that are perpendicular to the beam.
    # Gives a basis for refining beam
    axis = matrix.col((1, 0, 0))
    beamr0 = s0.cross(axis).normalize()
    beamr1 = beamr0.cross(s0).normalize()
    beamr2 = beamr1.cross(s0).normalize()

    assert approx_equal(s0.dot(beamr1), 0.0)
    assert approx_equal(s0.dot(beamr2), 0.0)
    assert approx_equal(beamr2.dot(beamr1), 0.0)
    # so the orthonormal vectors are s0, beamr1 and beamr2

    if mm_search_scope:
        plot_px_sz = experiments[0].detector[0].get_pixel_size()[0]
        plot_px_sz *= wide_search_binning
        grid = max(1, int(mm_search_scope / plot_px_sz))
        widegrid = 2 * grid + 1

        def get_experiment_score_for_coord(x, y):
            new_origin_offset = x * plot_px_sz * beamr1 + y * plot_px_sz * beamr2
            return sum(
                _get_origin_offset_score(
                    new_origin_offset,
                    solution_lists[i],
                    amax_lists[i],
                    reflection_lists[i],
                    experiment,
                )
                for i, experiment in enumerate(experiments)
            )

        scores = flex.double(
            get_experiment_score_for_coord(x, y)
            for y in range(-grid, grid + 1)
            for x in range(-grid, grid + 1)
        )

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

    # Do a simplex minimization
    class simplex_minimizer:
        def __init__(self, wide_search_offset):
            self.n = 2
            self.wide_search_offset = wide_search_offset
            self.optimizer = simplex_opt(
                dimension=self.n,
                matrix=[flex.random_double(self.n) for _ in range(self.n + 1)],
                evaluator=self,
                tolerance=1e-7,
            )
            self.x = self.optimizer.get_solution()
            self.offset = self.x[0] * 0.2 * beamr1 + self.x[1] * 0.2 * beamr2
            if self.wide_search_offset is not None:
                self.offset += self.wide_search_offset

        def target(self, vector):
            trial_origin_offset = vector[0] * 0.2 * beamr1 + vector[1] * 0.2 * beamr2
            if self.wide_search_offset is not None:
                trial_origin_offset += self.wide_search_offset
            target = 0
            for i, experiment in enumerate(experiments):
                target -= _get_origin_offset_score(
                    trial_origin_offset,
                    solution_lists[i],
                    amax_lists[i],
                    reflection_lists[i],
                    experiment,
                )
            return target

    new_offset = simplex_minimizer(wide_search_offset).offset

    if plot_search_scope:
        plot_px_sz = experiments[0].get_detector()[0].get_pixel_size()[0]
        grid = max(1, int(mm_search_scope / plot_px_sz))
        scores = flex.double()
        for y in range(-grid, grid + 1):
            for x in range(-grid, grid + 1):
                new_origin_offset = x * plot_px_sz * beamr1 + y * plot_px_sz * beamr2
                score = 0
                for i, experiment in enumerate(experiments):
                    score += _get_origin_offset_score(
                        new_origin_offset,
                        solution_lists[i],
                        amax_lists[i],
                        reflection_lists[i],
                        experiment,
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
            trial_origin_offset = (
                (idxs[idx_max % widegrid]) * beamr1
                + (idxs[idx_max // widegrid]) * beamr2
            )
            return trial_origin_offset

        show_plot(widegrid=2 * grid + 1, excursi=scores)

    new_experiments = copy.deepcopy(experiments)
    for expt in new_experiments:
        expt.detector = dps_extended.get_new_detector(expt.detector, new_offset)
    return new_experiments


def _get_origin_offset_score(
    trial_origin_offset, solutions, amax, spots_mm, experiment
):
    trial_detector = dps_extended.get_new_detector(
        experiment.detector, trial_origin_offset
    )
    experiment = copy.copy(experiment)
    experiment.detector = trial_detector

    # Key point for this is that the spots must correspond to detector
    # positions not to the correct RS position => reset any fixed rotation
    # to identity

    experiment.goniometer.set_fixed_rotation((1, 0, 0, 0, 1, 0, 0, 0, 1))
    spots_mm.map_centroids_to_reciprocal_space([experiment])
    return _sum_score_detail(spots_mm["rlp"], solutions, amax=amax)


def _sum_score_detail(reciprocal_space_vectors, solutions, granularity=None, amax=None):
    """Evaluates the probability that the trial value of
    (S0_vector | origin_offset) is correct, given the current estimate and the
    observations. The trial value comes through the reciprocal space vectors,
    and the current estimate comes through the short list of DPS solutions.
    Actual return value is a sum of NH terms, one for each DPS solution,
    each ranging from -1.0 to 1.0"""

    nh = min(solutions.size(), 20)  # extended API
    sum_score = 0.0
    for t in range(nh):
        dfft = Directional_FFT(
            angle=Direction(solutions[t]),
            xyzdata=reciprocal_space_vectors,
            granularity=5.0,
            amax=amax,  # extended API XXX
            F0_cutoff=11,
        )
        kval = dfft.kval()
        kmax = dfft.kmax()
        kval_cutoff = reciprocal_space_vectors.size() / 4.0
        if kval > kval_cutoff:
            ff = dfft.fft_result
            kbeam = ((-dfft.pmin) / dfft.delta_p) + 0.5
            Tkmax = cmath.phase(ff[kmax])
            backmax = math.cos(
                Tkmax + (2 * math.pi * kmax * kbeam / (2 * ff.size() - 1))
            )
            # Here it should be possible to calculate a gradient.
            # Then minimize with respect to two coordinates.
            # Use lbfgs?  Have second derivatives?
            # can I do something local to model the cosine wave?
            # direction of wave travel.  Period. phase.
            sum_score += backmax
    return sum_score


def run_dps(experiment, spots_mm, max_cell):
    # max_cell: max possible cell in Angstroms; set to None,
    # determine from data
    # recommended_grid_sampling_rad: grid sampling in radians; guess for now

    horizon_phil = iotbx.phil.parse(input_string=indexing_api_defs).extract()
    DPS = DPS_primitive_lattice(
        max_cell=max_cell, recommended_grid_sampling_rad=None, horizon_phil=horizon_phil
    )

    DPS.S0_vector = matrix.col(experiment.beam.get_s0())
    DPS.inv_wave = 1.0 / experiment.beam.get_wavelength()
    if experiment.goniometer is None:
        DPS.axis = matrix.col((1, 0, 0))
    else:
        DPS.axis = matrix.col(experiment.goniometer.get_rotation_axis())
    DPS.set_detector(experiment.detector)

    # transform input into what DPS needs
    # i.e., construct a flex.vec3 double consisting of mm spots, phi in degrees
    data = flex.vec3_double()
    for spot in spots_mm.rows():
        data.append(
            (
                spot["xyzobs.mm.value"][0],
                spot["xyzobs.mm.value"][1],
                spot["xyzobs.mm.value"][2] * 180.0 / math.pi,
            )
        )

    logger.info("Running DPS using %i reflections", len(data))

    DPS.index(
        raw_spot_input=data,
        panel_addresses=flex.int(s["panel"] for s in spots_mm.rows()),
    )
    solutions = DPS.getSolutions()

    logger.info(
        "Found %i solution%s with max unit cell %.2f Angstroms.",
        len(solutions),
        plural_s(len(solutions))[1],
        DPS.amax,
    )

    # There must be at least 3 solutions to make a set,
    # otherwise return empty result
    if len(solutions) < 3:
        return {}
    return {"solutions": flex.vec3_double(s.dvec for s in solutions), "amax": DPS.amax}


def discover_better_experimental_model(
    experiments,
    reflections,
    params,
    nproc=1,
    d_min=None,
    mm_search_scope=4.0,
    wide_search_binning=1,
    plot_search_scope=False,
):
    assert len(experiments) == len(reflections)
    assert len(experiments) > 0

    refl_lists = []
    max_cell_list = []

    # The detector/beam of the first experiment is used to define the basis
    # for the optimisation, so assert that the beam intersects with
    # the detector
    detector = experiments[0].detector
    beam = experiments[0].beam
    beam_panel = detector.get_panel_intersection(beam.get_s0())
    if beam_panel == -1:
        raise Sorry("input beam does not intersect detector")

    for expt, refl in zip(experiments, reflections):
        refl = copy.deepcopy(refl)
        refl["imageset_id"] = flex.int(len(refl), 0)
        refl.centroid_px_to_mm([expt])
        refl.map_centroids_to_reciprocal_space([expt])

        if d_min is not None:
            d_spacings = 1 / refl["rlp"].norms()
            sel = d_spacings > d_min
            refl = refl.select(sel)

        # derive a max_cell from mm spots
        if params.max_cell is None:
            max_cell = find_max_cell(
                refl, max_cell_multiplier=1.3, step_size=45
            ).max_cell
            max_cell_list.append(max_cell)

        refl_condition = refl.size() > params.max_reflections
        if params.max_reflections is not None and refl_condition:
            logger.info(
                "Selecting subset of %i reflections for analysis",
                params.max_reflections,
            )
            perm = flex.random_permutation(refl.size())
            sel = perm[: params.max_reflections]
            refl = refl.select(sel)

        refl_lists.append(refl)

    if params.max_cell is None:
        max_cell = flex.median(flex.double(max_cell_list))
    else:
        max_cell = params.max_cell

    with concurrent.futures.ProcessPoolExecutor(max_workers=nproc) as pool:
        solution_lists = []
        amax_list = []
        for result in pool.map(
            run_dps, experiments, refl_lists, itertools.repeat(max_cell)
        ):
            if result.get("solutions"):
                solution_lists.append(result["solutions"])
                amax_list.append(result["amax"])

    if not solution_lists:
        raise Sorry("No solutions found")

    new_experiments = optimize_origin_offset_local_scope(
        experiments,
        refl_lists,
        solution_lists,
        amax_list,
        mm_search_scope=mm_search_scope,
        wide_search_binning=wide_search_binning,
        plot_search_scope=plot_search_scope,
    )
    new_detector = new_experiments[0].detector
    old_panel, old_beam_centre = detector.get_ray_intersection(beam.get_s0())
    (new_panel, new_beam_centre) = new_detector.get_ray_intersection(beam.get_s0())

    old_bc_px = detector[old_panel].millimeter_to_pixel(old_beam_centre)
    new_bc_px = new_detector[new_panel].millimeter_to_pixel(new_beam_centre)

    logger.info(
        "Old beam centre: {:.2f}, {:.2f} mm".format(*old_beam_centre)
        + " ({:.1f}, {:.1f} px)".format(*old_bc_px)
    )
    logger.info(
        "New beam centre: {:.2f}, {:.2f} mm".format(*new_beam_centre)
        + " ({:.1f}, {:.1f} px)".format(*new_bc_px)
    )
    logger.info(
        "Shift: {:.2f}, {:.2f} mm".format(
            *(matrix.col(old_beam_centre) - matrix.col(new_beam_centre)).elems
        )
        + " ({:.1f}, {:.1f} px)".format(
            *(matrix.col(old_bc_px) - matrix.col(new_bc_px)).elems
        )
    )
    return new_experiments
