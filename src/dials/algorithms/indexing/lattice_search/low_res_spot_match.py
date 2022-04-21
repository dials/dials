from __future__ import annotations

import copy
import logging
import math
import operator

import libtbx.phil
from cctbx import miller
from dxtbx.model import Crystal
from scitbx import matrix
from scitbx.math import least_squares_plane, superpose

from dials.algorithms.indexing import DialsIndexError
from dials.array_family import flex

from .strategy import Strategy

logger = logging.getLogger(__name__)

TWO_PI = 2.0 * math.pi
FIVE_DEG = TWO_PI * 5.0 / 360.0


class CompleteGraph:
    def __init__(self, seed_vertex):
        self.vertices = [seed_vertex]
        self.weight = [{0: 0.0}]
        self.total_weight = 0.0

    def factory_add_vertex(self, vertex, weights_to_other):
        # Return a new graph as a copy of this with an extra vertex added. This
        # is a factory rather than a change in-place because CompleteGraph ought
        # to be immutable to implement __hash__
        g = copy.deepcopy(self)

        current_len = len(g.vertices)
        assert len(weights_to_other) == current_len
        g.vertices.append(vertex)
        node = current_len

        # Update distances from other nodes to the new one
        for i, w in enumerate(weights_to_other):
            g.weight[i][node] = w

        # Add distances to other nodes from this one
        weights_to_other.append(0.0)
        to_other = {}
        for i, w in enumerate(weights_to_other):
            to_other[i] = w
        g.weight.append(to_other)

        # Update the total weight
        g.total_weight += sum(weights_to_other)

        # Sort the vertices and weights by spot_id
        l = zip(g.vertices, g.weight)
        l = sorted(l, key=lambda v_w: v_w[0]["spot_id"])
        v, w = zip(*l)
        g.vertices = list(v)
        g.weight = list(w)

        return g

    def __hash__(self):
        h = tuple((e["spot_id"], e["miller_index"]) for e in self.vertices)
        return hash(h)

    def __eq__(self, other):
        for a, b in zip(self.vertices, other.vertices):
            if a["spot_id"] != b["spot_id"]:
                return False
            if a["miller_index"] != b["miller_index"]:
                return False
        return True

    def __ne__(self, other):
        return not self == other


low_res_spot_match_phil_str = """\
candidate_spots
{
    limit_resolution_by = *n_spots d_min
    .type = choice

    d_min = 15.0
    .type = float(value_min=0)

    n_spots = 10
    .type = int

    d_star_tolerance = 4.0
    .help = "Number of sigmas from the centroid position for which to "
            "calculate d* bands"
    .type = float
}

use_P1_indices_as_seeds = False
    .type = bool

search_depth = *triplets quads
    .type = choice

bootstrap_crystal = False
    .type = bool

max_pairs = 200
    .type = int

max_triplets = 600
    .type = int

max_quads = 600
    .type = int
"""


class LowResSpotMatch(Strategy):
    """
    A lattice search strategy matching low resolution spots to candidate indices.

    The match is based on resolution and reciprocal space distance between observed
    spots. A prior unit cell and space group are required and solutions are assessed
    by matching the low resolution spots against candidate reflection positions
    predicted from the known cell. This lattice search strategy is a special case
    designed to work for electron diffraction still images, in which one typically
    only collects reflections from the zero-order Laue zone. In principle, it is not
    limited to this type of data, but probably works best with narrow wedges,
    good initial geometry and a small beam-stop shadow so that a good number of
    low-order reflections are collected.
    """

    phil_help = (
        "A lattice search strategy that matches low resolution spots to candidate "
        "indices based on a known unit cell and space group. Designed primarily for "
        "electron diffraction still images."
    )

    phil_scope = libtbx.phil.parse(low_res_spot_match_phil_str)

    def __init__(
        self, target_symmetry_primitive, max_lattices, params=None, *args, **kwargs
    ):
        """Construct a LowResSpotMatch object.

        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target
                crystal symmetry and unit cell

            max_lattices (int): The maximum number of lattice models to find
        """
        super().__init__(params=params, *args, **kwargs)
        self._target_symmetry_primitive = target_symmetry_primitive
        self._max_lattices = max_lattices

        if target_symmetry_primitive is None:
            raise DialsIndexError(
                "Target unit cell and space group must be provided for low_res_spot_match"
            )

        # Set reciprocal space orthogonalisation matrix
        uc = self._target_symmetry_primitive.unit_cell()
        self.Bmat = matrix.sqr(uc.fractionalization_matrix()).transpose()

    def find_crystal_models(self, reflections, experiments):
        """Find a list of candidate crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models
        """

        # Take a subset of the observations at the same resolution and calculate
        # some values that will be needed for the search
        self._calc_obs_data(reflections, experiments)

        # Construct a library of candidate low res indices with their d* values
        self._calc_candidate_hkls()

        # First search: match each observation with candidate indices within the
        # acceptable resolution band
        self._calc_seeds_and_stems()
        if self._params.use_P1_indices_as_seeds:
            seeds = self.stems
        else:
            seeds = self.seeds
        logger.info("Using %s seeds", len(seeds))

        # Second search: match seed spots with another spot from a different
        # reciprocal lattice row, such that the observed reciprocal space distances
        # are within tolerances
        pairs = []
        for seed in seeds:
            pairs.extend(self._pairs_with_seed(seed))
        logger.info("Found %s pairs", len(pairs))
        pairs = list(set(pairs))  # filter duplicates

        if self._params.max_pairs:
            pairs.sort(key=operator.attrgetter("total_weight"))
            idx = self._params.max_pairs
            pairs = pairs[0:idx]
        logger.info("Using %s highest-scoring pairs", len(pairs))

        # Further search iterations: extend to more spots within tolerated distances
        triplets = []
        for pair in pairs:
            triplets.extend(self._extend_by_candidates(pair))
        logger.info("Found %s triplets", len(triplets))
        triplets = list(set(triplets))  # filter duplicates
        if self._params.max_triplets:
            triplets.sort(key=operator.attrgetter("total_weight"))
            idx = self._params.max_triplets
            triplets = triplets[0:idx]
        logger.info("Using %s highest-scoring triplets", len(triplets))

        branches = triplets
        if self._params.search_depth == "quads":
            quads = []
            for triplet in triplets:
                quads.extend(self._extend_by_candidates(triplet))
            logger.info("%s quads", len(quads))
            quads = list(set(quads))  # filter duplicates
            if self._params.max_quads:
                quads.sort(key=operator.attrgetter("total_weight"))
                idx = self._params.max_quads
                quads = quads[0:idx]
            logger.info("Using %s highest-scoring quads", len(quads))
            branches = quads

        # Sort branches by total deviation of observed distances from expected
        branches.sort(key=operator.attrgetter("total_weight"))

        candidate_crystal_models = []
        for branch in branches:
            model = self._fit_crystal_model(branch)
            if model:
                candidate_crystal_models.append(model)
            if len(candidate_crystal_models) == self._max_lattices:
                break

        self.candidate_crystal_models = candidate_crystal_models
        return self.candidate_crystal_models

    def _calc_candidate_hkls(self):
        # First a list of indices that fill 1 ASU
        hkl_list = miller.build_set(
            self._target_symmetry_primitive,
            anomalous_flag=False,
            d_min=self._params.candidate_spots.d_min,
        )
        rt = flex.reflection_table()
        rt["miller_index"] = hkl_list.indices()
        rt["d_star"] = 1.0 / hkl_list.d_spacings().data()
        rt["rlp_datum"] = self.Bmat.elems * rt["miller_index"].as_vec3_double()
        self.candidate_hkls = rt

        # Now P1 indices with separate Friedel pairs
        hkl_list = miller.build_set(
            self._target_symmetry_primitive,
            anomalous_flag=True,
            d_min=self._params.candidate_spots.d_min,
        )
        hkl_list_p1 = hkl_list.expand_to_p1()
        rt = flex.reflection_table()
        rt["miller_index"] = hkl_list_p1.indices()
        rt["d_star"] = 1.0 / hkl_list_p1.d_spacings().data()
        rt["rlp_datum"] = self.Bmat.elems * rt["miller_index"].as_vec3_double()
        self.candidate_hkls_p1 = rt
        return

    def _calc_obs_data(self, reflections, experiments):
        """Calculates a set of low resolution observations to try to match to
        indices. Each observation will record its d* value as well as
        tolerated d* bands and a 'clock angle'"""

        spot_d_star = reflections["rlp"].norms()
        if self._params.candidate_spots.limit_resolution_by == "n_spots":
            n_spots = self._params.candidate_spots.n_spots
            n_spots = min(n_spots, len(reflections) - 1)
            d_star_max = flex.sorted(spot_d_star)[n_spots - 1]
            self._params.candidate_spots.d_min = 1.0 / d_star_max

        # First select low resolution spots only
        spot_d_star = reflections["rlp"].norms()
        d_star_max = 1.0 / self._params.candidate_spots.d_min
        sel = spot_d_star <= d_star_max
        self.spots = reflections.select(sel)
        self.spots["d_star"] = spot_d_star.select(sel)

        # XXX In what circumstance might there be more than one experiment?
        detector = experiments.detectors()[0]
        beam = experiments.beams()[0]

        # Lab coordinate of the beam centre, using the first spot's panel
        panel = detector[self.spots[0]["panel"]]
        bc = panel.get_ray_intersection(beam.get_s0())
        bc_lab = panel.get_lab_coord(bc)

        # Lab coordinate of each spot
        spot_lab = flex.vec3_double(len(self.spots))
        pnl_ids = set(self.spots["panel"])
        for pnl in pnl_ids:
            sel = self.spots["panel"] == pnl
            panel = detector[pnl]
            obs = self.spots["xyzobs.mm.value"].select(sel)
            x_mm, y_mm, _ = obs.parts()
            spot_lab.set_selected(
                sel, panel.get_lab_coord(flex.vec2_double(x_mm, y_mm))
            )

        # Radius vectors for each spot
        radius = spot_lab - bc_lab

        # Usually the radius vectors would all be in a single plane, but this might
        # not be the case if the spots are on different panels. To put them on the
        # same plane, project onto fast/slow of the panel used to get the beam
        # centre
        df = flex.vec3_double(len(self.spots), detector[0].get_fast_axis())
        ds = flex.vec3_double(len(self.spots), detector[0].get_slow_axis())
        clock_dirs = (radius.dot(df) * df + radius.dot(ds) * ds).each_normalize()

        # From this, find positive angles of each vector around a clock, using the
        # fast axis as 12 o'clock
        angs = clock_dirs.angle(detector[0].get_fast_axis())
        dots = clock_dirs.dot(detector[0].get_slow_axis())
        sel = dots < 0  # select directions in the second half of the clock face
        angs.set_selected(sel, (TWO_PI - angs.select(sel)))
        self.spots["clock_angle"] = angs

        # Project radius vectors onto fast/slow of the relevant panels
        df = flex.vec3_double(len(self.spots))
        ds = flex.vec3_double(len(self.spots))
        for pnl in pnl_ids:
            sel = self.spots["panel"] == pnl
            panel = detector[pnl]
            df.set_selected(sel, panel.get_fast_axis())
            ds.set_selected(sel, panel.get_slow_axis())
        panel_dirs = (radius.dot(df) * df + radius.dot(ds) * ds).each_normalize()

        # Calc error along each panel direction with simple error propagation
        # that assumes no covariance between x and y centroid errors.
        x = panel_dirs.dot(df)
        y = panel_dirs.dot(ds)
        x2, y2 = flex.pow2(x), flex.pow2(y)
        r2 = x2 + y2
        sig_x2, sig_y2, _ = self.spots["xyzobs.mm.variance"].parts()
        var_r = (x2 / r2) * sig_x2 + (y2 / r2) * sig_y2
        sig_r = flex.sqrt(var_r)

        # Pixel coordinates at limits of the band
        tol = self._params.candidate_spots.d_star_tolerance
        outer_spot_lab = spot_lab + panel_dirs * (tol * sig_r)
        inner_spot_lab = spot_lab - panel_dirs * (tol * sig_r)

        # Set d* at band limits
        inv_lambda = 1.0 / beam.get_wavelength()
        s1_outer = outer_spot_lab.each_normalize() * inv_lambda
        s1_inner = inner_spot_lab.each_normalize() * inv_lambda
        self.spots["d_star_outer"] = (s1_outer - beam.get_s0()).norms()
        self.spots["d_star_inner"] = (s1_inner - beam.get_s0()).norms()
        self.spots["d_star_band2"] = flex.pow2(
            self.spots["d_star_outer"] - self.spots["d_star_inner"]
        )

    def _calc_seeds_and_stems(self):
        # As the first stage of search, determine a list of seed spots for further
        # stages. Order these by distance of observed d* from the candidate
        # reflection's canonical d*

        # First the 'seeds' (in 1 ASU)
        self.seeds = []
        for i, spot in enumerate(self.spots.rows()):
            sel = (self.candidate_hkls["d_star"] <= spot["d_star_outer"]) & (
                self.candidate_hkls["d_star"] >= spot["d_star_inner"]
            )
            cands = self.candidate_hkls.select(sel)
            for c in cands.rows():
                r_dst = abs(c["d_star"] - spot["d_star"])
                self.seeds.append(
                    {
                        "spot_id": i,
                        "miller_index": c["miller_index"],
                        "rlp_datum": matrix.col(c["rlp_datum"]),
                        "residual_d_star": r_dst,
                        "clock_angle": spot["clock_angle"],
                    }
                )

        self.seeds.sort(key=operator.itemgetter("residual_d_star"))

        # Now the 'stems' to use in second search level, using all indices in P 1
        self.stems = []
        for i, spot in enumerate(self.spots.rows()):
            sel = (self.candidate_hkls_p1["d_star"] <= spot["d_star_outer"]) & (
                self.candidate_hkls_p1["d_star"] >= spot["d_star_inner"]
            )
            cands = self.candidate_hkls_p1.select(sel)
            for c in cands.rows():
                r_dst = abs(c["d_star"] - spot["d_star"])
                self.stems.append(
                    {
                        "spot_id": i,
                        "miller_index": c["miller_index"],
                        "rlp_datum": matrix.col(c["rlp_datum"]),
                        "residual_d_star": r_dst,
                        "clock_angle": spot["clock_angle"],
                    }
                )

        self.stems.sort(key=operator.itemgetter("residual_d_star"))

    def _pairs_with_seed(self, seed):
        seed_rlp = matrix.col(self.spots[seed["spot_id"]]["rlp"])

        result = []
        for cand in self.stems:
            # Don't check the seed spot itself
            if cand["spot_id"] == seed["spot_id"]:
                continue

            # Skip spots at a very similar clock angle, which probably belong to the
            # same line of indices from the origin
            angle_diff = cand["clock_angle"] - seed["clock_angle"]
            angle_diff = abs(((angle_diff + math.pi) % TWO_PI) - math.pi)
            if angle_diff < FIVE_DEG:
                continue

            # Calculate the plane normal for the plane containing the seed and stem.
            # Skip pairs of Miller indices that belong to the same line
            seed_vec = seed["rlp_datum"]
            cand_vec = cand["rlp_datum"]
            try:
                seed_vec.cross(cand_vec).normalize()
            except ZeroDivisionError:
                continue

            # Compare expected reciprocal space distance with observed distance
            cand_rlp = matrix.col(self.spots[cand["spot_id"]]["rlp"])
            obs_dist = (cand_rlp - seed_rlp).length()
            exp_dist = (seed_vec - cand_vec).length()
            r_dist = abs(obs_dist - exp_dist)

            # If the distance difference is larger than the sum in quadrature of the
            # tolerated d* bands then reject the candidate
            sq_band1 = self.spots[seed["spot_id"]]["d_star_band2"]
            sq_band2 = self.spots[cand["spot_id"]]["d_star_band2"]
            if r_dist > math.sqrt(sq_band1 + sq_band2):
                continue

            # Store the seed-stem match as a 2-node graph
            g = CompleteGraph(
                {
                    "spot_id": seed["spot_id"],
                    "miller_index": seed["miller_index"],
                    "rlp_datum": seed["rlp_datum"],
                }
            )
            g = g.factory_add_vertex(
                {
                    "spot_id": cand["spot_id"],
                    "miller_index": cand["miller_index"],
                    "rlp_datum": cand["rlp_datum"],
                },
                weights_to_other=[r_dist],
            )
            result.append(g)
        return result

    def _extend_by_candidates(self, graph):

        existing_ids = [e["spot_id"] for e in graph.vertices]
        obs_relps = [matrix.col(self.spots[e]["rlp"]) for e in existing_ids]
        exp_relps = [e["rlp_datum"] for e in graph.vertices]

        result = []

        for cand in self.stems:
            # Don't check spots already matched
            if cand["spot_id"] in existing_ids:
                continue

            # Compare expected reciprocal space distances with observed distances
            cand_rlp = matrix.col(self.spots[cand["spot_id"]]["rlp"])
            cand_vec = cand["rlp_datum"]

            obs_dists = [(cand_rlp - rlp).length() for rlp in obs_relps]
            exp_dists = [(vec - cand_vec).length() for vec in exp_relps]

            residual_dist = [abs(a - b) for (a, b) in zip(obs_dists, exp_dists)]

            # If any of the distance differences is larger than the sum in quadrature
            # of the tolerated d* bands then reject the candidate
            sq_candidate_band = self.spots[cand["spot_id"]]["d_star_band2"]
            bad_candidate = False
            for r_dist, spot_id in zip(residual_dist, existing_ids):
                sq_relp_band = self.spots[spot_id]["d_star_band2"]
                if r_dist > math.sqrt(sq_relp_band + sq_candidate_band):
                    bad_candidate = True
                    break
            if bad_candidate:
                continue

            # Calculate co-planarity of the relps, including the origin
            points = flex.vec3_double(exp_relps + [cand_vec, (0.0, 0.0, 0.0)])
            plane = least_squares_plane(points)
            plane_score = flex.sum_sq(
                points.dot(plane.normal) - plane.distance_to_origin
            )

            # Reject if the group of relps are too far from lying in a single plane.
            # This cut-off was determined by trial and error using simulated images.
            if plane_score > 6e-7:
                continue

            # Construct a graph including the accepted candidate node
            g = graph.factory_add_vertex(
                {
                    "spot_id": cand["spot_id"],
                    "miller_index": cand["miller_index"],
                    "rlp_datum": cand["rlp_datum"],
                },
                weights_to_other=residual_dist,
            )

            result.append(g)

        return result

    @staticmethod
    def _fit_U_from_superposed_points(reference, other):

        # Add the origin to both sets of points
        reference.append((0, 0, 0))
        other.append((0, 0, 0))

        # Find U matrix that takes ideal relps to the reference
        fit = superpose.least_squares_fit(reference, other)
        return fit.r

    def _fit_crystal_model(self, graph):

        vertices = graph.vertices

        # Reciprocal lattice points of the observations
        sel = flex.size_t([e["spot_id"] for e in vertices])
        reference = self.spots["rlp"].select(sel)

        # Ideal relps from the known cell
        other = flex.vec3_double([e["rlp_datum"] for e in vertices])

        U = self._fit_U_from_superposed_points(reference, other)
        UB = U * self.Bmat

        if self._params.bootstrap_crystal:

            # Attempt to index the low resolution spots
            from dials_algorithms_indexing_ext import AssignIndices

            phi = self.spots["xyzobs.mm.value"].parts()[2]
            UB_matrices = flex.mat3_double([UB])
            result = AssignIndices(self.spots["rlp"], phi, UB_matrices, tolerance=0.3)
            hkl = result.miller_indices()
            sel = hkl != (0, 0, 0)
            hkl_vec = hkl.as_vec3_double().select(sel)

            # Use the result to get a new UB matrix
            reference = self.spots["rlp"].select(sel)
            other = self.Bmat.elems * hkl_vec
            U = self._fit_U_from_superposed_points(reference, other)
            UB = U * self.Bmat

        # Calculate RMSD of the fit
        rms = reference.rms_difference(U.elems * other)

        # Construct a crystal model
        xl = Crystal(A=UB, space_group_symbol="P1")

        # Monkey-patch crystal to return rms of the fit (useful?)
        xl.rms = rms

        return xl
