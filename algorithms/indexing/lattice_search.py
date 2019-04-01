#!/usr/bin/env python
# -*- mode: python; coding: utf-8; indent-tabs-mode: nil; python-indent: 2 -*-
#
# dials.algorithms.indexing.fft1d.py
#
#  Copyright (C) 2014 Diamond Light Source
#
#  Author: Richard Gildea
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
from __future__ import print_function

import logging
import math

from scitbx.array_family import flex
import scitbx.matrix
from dxtbx.model.experiment_list import Experiment, ExperimentList

from dials.algorithms.indexing.basis_vector_search import strategies
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)
from dials.algorithms.indexing.basis_vector_search import optimise

logger = logging.getLogger(__name__)


from dials.algorithms.indexing import indexer


class BasisVectorSearch(indexer.indexer_base):
    def __init__(self, reflections, experiments, params=None):
        super(BasisVectorSearch, self).__init__(reflections, experiments, params)
        if params.indexing.method == "fft1d":
            self._basis_vector_search_strategy = strategies.fft1d(
                max_cell=self.params.max_cell,
                characteristic_grid=self.params.fft1d.characteristic_grid,
            )
        elif params.indexing.method == "fft3d":
            self._basis_vector_search_strategy = strategies.fft3d(
                self.params.max_cell,
                self.params.fft3d.reciprocal_space_grid.n_points,
                d_min=self.params.fft3d.reciprocal_space_grid.d_min,
                b_iso=self.params.b_iso,
                rmsd_cutoff=self.params.rmsd_cutoff,
                peak_volume_cutoff=self.params.fft3d.peak_volume_cutoff,
                min_cell=self.params.min_cell,
            )
        elif params.indexing.method == "real_space_grid_search":
            assert self._symmetry_handler.target_symmetry_primitive is not None
            assert (
                self._symmetry_handler.target_symmetry_primitive.unit_cell() is not None
            )
            self._basis_vector_search_strategy = strategies.real_space_grid_search(
                max_cell=self.params.max_cell,
                target_unit_cell=self._symmetry_handler.target_symmetry_primitive.unit_cell(),
                characteristic_grid=self.params.real_space_grid_search.characteristic_grid,
            )
        else:
            raise RuntimeError("Unknown basis vector search method")

    def find_candidate_basis_vectors(self):
        self.d_min = self.params.refinement_protocol.d_min_start
        sel = self.reflections["id"] == -1
        if self.d_min is not None:
            sel &= 1 / self.reflections["rlp"].norms() > self.d_min
        reflections = self.reflections.select(sel)
        self.candidate_basis_vectors, used_in_indexing = self._basis_vector_search_strategy.find_basis_vectors(
            reflections["rlp"]
        )
        self._used_in_indexing = sel.iselection().select(used_in_indexing)
        if self.d_min is None:
            self.d_min = flex.min(
                1 / self.reflections["rlp"].select(self._used_in_indexing).norms()
            )

        self.debug_show_candidate_basis_vectors()
        if self.params.debug_plots:
            self.debug_plot_candidate_basis_vectors()

        return self.candidate_basis_vectors

    def find_lattices(self):
        self.find_candidate_basis_vectors()
        if self.params.optimise_initial_basis_vectors:
            self.optimise_basis_vectors()
        self.candidate_crystal_models = self.find_candidate_orientation_matrices(
            self.candidate_basis_vectors
        )
        crystal_model, n_indexed = self.choose_best_orientation_matrix(
            self.candidate_crystal_models
        )
        if crystal_model is not None:
            crystal_models = [crystal_model]
        else:
            crystal_models = []
        experiments = ExperimentList()
        for cm in crystal_models:
            for expt in self.experiments:
                experiments.append(
                    Experiment(
                        imageset=expt.imageset,
                        beam=expt.beam,
                        detector=expt.detector,
                        goniometer=expt.goniometer,
                        scan=expt.scan,
                        crystal=cm,
                    )
                )
        return experiments

    def find_candidate_orientation_matrices(self, candidate_basis_vectors):
        from dials.algorithms.indexing.basis_vector_search import combinations

        candidate_crystal_models = combinations.candidate_orientation_matrices(
            candidate_basis_vectors
        )
        if self._symmetry_handler.target_symmetry_reference_setting is not None:
            target_symmetry = self._symmetry_handler.target_symmetry_reference_setting
        elif self._symmetry_handler.target_symmetry_primitive is not None:
            target_symmetry = self._symmetry_handler.target_symmetry_primitive
        else:
            target_symmetry = None
        if target_symmetry is not None:
            candidate_crystal_models = combinations.filter_known_symmetry(
                candidate_crystal_models,
                target_symmetry,
                relative_length_tolerance=self.params.known_symmetry.relative_length_tolerance,
                absolute_angle_tolerance=self.params.known_symmetry.absolute_angle_tolerance,
                max_delta=self.params.known_symmetry.max_delta,
            )

        candidate_crystal_models = self.filter_similar_orientations(
            candidate_crystal_models
        )

        return candidate_crystal_models

    def choose_best_orientation_matrix(self, candidate_orientation_matrices):

        from dials.algorithms.indexing import model_evaluation

        solution_scorer = self.params.basis_vector_combinations.solution_scorer
        if solution_scorer == "weighted":
            weighted_params = self.params.basis_vector_combinations.weighted
            solutions = model_evaluation.ModelRankWeighted(
                power=weighted_params.power,
                volume_weight=weighted_params.volume_weight,
                n_indexed_weight=weighted_params.n_indexed_weight,
                rmsd_weight=weighted_params.rmsd_weight,
            )
        else:
            filter_params = self.params.basis_vector_combinations.filter
            solutions = model_evaluation.ModelRankFilter(
                check_doubled_cell=filter_params.check_doubled_cell,
                likelihood_cutoff=filter_params.likelihood_cutoff,
                volume_cutoff=filter_params.volume_cutoff,
                n_indexed_cutoff=filter_params.n_indexed_cutoff,
            )

        args = []

        for cm in candidate_orientation_matrices:
            sel = self.reflections["id"] == -1
            if self.d_min is not None:
                sel &= 1 / self.reflections["rlp"].norms() > self.d_min
            xo, yo, zo = self.reflections["xyzobs.mm.value"].parts()
            imageset_id = self.reflections["imageset_id"]
            experiments = ExperimentList()
            for i_expt, expt in enumerate(self.experiments):
                # XXX Not sure if we still need this loop over self.experiments
                if expt.scan is not None:
                    start, end = expt.scan.get_oscillation_range()
                    if (end - start) > 360:
                        # only use reflections from the first 360 degrees of the scan
                        sel.set_selected(
                            (imageset_id == i_expt)
                            & (zo > ((start * math.pi / 180) + 2 * math.pi)),
                            False,
                        )
                experiments.append(
                    Experiment(
                        imageset=expt.imageset,
                        beam=expt.beam,
                        detector=expt.detector,
                        goniometer=expt.goniometer,
                        scan=expt.scan,
                        crystal=cm,
                    )
                )
            refl = self.reflections.select(sel)
            self.index_reflections(experiments, refl)
            if refl.get_flags(refl.flags.indexed).count(True) == 0:
                continue

            from rstbx.dps_core.cell_assessment import SmallUnitCellVolume
            from dials.algorithms.indexing import non_primitive_basis

            threshold = self.params.basis_vector_combinations.sys_absent_threshold
            if threshold and (
                self._symmetry_handler.target_symmetry_primitive is None
                or self._symmetry_handler.target_symmetry_primitive.unit_cell() is None
            ):
                try:
                    non_primitive_basis.correct(
                        experiments, refl, self._assign_indices, threshold
                    )
                    if refl.get_flags(refl.flags.indexed).count(True) == 0:
                        continue
                except SmallUnitCellVolume:
                    logger.debug(
                        "correct_non_primitive_basis SmallUnitCellVolume error for unit cell %s:"
                        % experiments[0].crystal.get_unit_cell()
                    )
                    continue
                except RuntimeError as e:
                    if "Krivy-Gruber iteration limit exceeded" in str(e):
                        logger.debug(
                            "correct_non_primitive_basis Krivy-Gruber iteration limit exceeded error for unit cell %s:"
                            % experiments[0].crystal.get_unit_cell()
                        )
                        continue
                    raise
                if (
                    experiments[0].crystal.get_unit_cell().volume()
                    < self.params.min_cell_volume
                ):
                    continue

            if self.params.known_symmetry.space_group is not None:
                new_crystal, cb_op_to_primitive = self._symmetry_handler.apply_symmetry(
                    experiments[0].crystal
                )
                if new_crystal is None:
                    continue
                experiments[0].crystal.update(new_crystal)
                if not cb_op_to_primitive.is_identity_op():
                    sel = refl["id"] > -1
                    miller_indices = refl["miller_index"].select(sel)
                    miller_indices = cb_op_to_primitive.apply(miller_indices)
                    refl["miller_index"].set_selected(sel, miller_indices)
                if 0 and self.cb_op_primitive_to_given is not None:
                    sel = refl["id"] > -1
                    experiments[0].crystal.update(
                        experiments[0].crystal.change_basis(
                            self.cb_op_primitive_to_given
                        )
                    )
                    miller_indices = refl["miller_index"].select(sel)
                    miller_indices = self.cb_op_primitive_to_given.apply(miller_indices)
                    refl["miller_index"].set_selected(sel, miller_indices)

            args.append((experiments, refl))
            if len(args) == self.params.basis_vector_combinations.max_refine:
                break

        from libtbx import easy_mp

        evaluator = model_evaluation.ModelEvaluation(self.all_params)
        results = easy_mp.parallel_map(
            evaluator.evaluate,
            args,
            iterable_type=easy_mp.posiargs,
            processes=self.params.nproc,
            preserve_exception_message=True,
        )

        for soln in results:
            if soln is None:
                continue
            solutions.append(soln)

        if len(solutions):
            logger.info("Candidate solutions:")
            logger.info(str(solutions))
            best_model = solutions.best_model()
            logger.debug("best model_likelihood: %.2f" % best_model.model_likelihood)
            logger.debug("best n_indexed: %i" % best_model.n_indexed)
            self.hkl_offset = best_model.hkl_offset
            return best_model.crystal, best_model.n_indexed
        else:
            return None, None

    def filter_similar_orientations(self, crystal_models):

        for cryst in crystal_models:

            if (
                self.refined_experiments is not None
                and len(self.refined_experiments) > 0
            ):

                orientation_too_similar = False
                for i_a, cryst_a in enumerate(self.refined_experiments.crystals()):
                    R_ab, axis, angle, cb_op_ab = difference_rotation_matrix_axis_angle(
                        cryst_a, cryst
                    )
                    min_angle = (
                        self.params.multiple_lattice_search.minimum_angular_separation
                    )
                    if abs(angle) < min_angle:  # degrees
                        orientation_too_similar = True
                        break
                if orientation_too_similar:
                    logger.debug("skipping crystal: too similar to other crystals")
                    continue
            yield cryst

    def optimise_basis_vectors(self):

        optimised_basis_vectors = optimise.optimise_basis_vectors(
            self.reflections["rlp"].select(self._used_in_indexing),
            self.candidate_basis_vectors,
        )
        self.candidate_basis_vectors = [
            scitbx.matrix.col(v) for v in optimised_basis_vectors
        ]

    def debug_show_candidate_basis_vectors(self):

        vectors = self.candidate_basis_vectors

        logger.debug("Candidate basis vectors:")
        for i, v in enumerate(vectors):
            logger.debug("%s %s" % (i, v.length()))  # , vector_heights[i]

        if self.params.debug:
            # print a table of the angles between each pair of vectors

            from cStringIO import StringIO

            s = StringIO()

            angles = flex.double(len(vectors) ** 2)
            angles.reshape(flex.grid(len(vectors), len(vectors)))

            for i in range(len(vectors)):
                v_i = vectors[i]
                for j in range(i + 1, len(vectors)):
                    v_j = vectors[j]
                    angles[i, j] = v_i.angle(v_j, deg=True)

            print((" " * 7), end=" ", file=s)
            for i in range(len(vectors)):
                print("%7.3f" % vectors[i].length(), end=" ", file=s)
            print(file=s)
            for i in range(len(vectors)):
                print("%7.3f" % vectors[i].length(), end=" ", file=s)
                for j in range(len(vectors)):
                    if j <= i:
                        print((" " * 7), end=" ", file=s)
                    else:
                        print("%5.1f  " % angles[i, j], end=" ", file=s)
                print(file=s)

            logger.debug(s.getvalue())

    def debug_plot_candidate_basis_vectors(self):
        from matplotlib import pyplot
        from mpl_toolkits.mplot3d import Axes3D  # noqa: F401

        fig = pyplot.figure()
        ax = fig.add_subplot(111, projection="3d")
        ax.scatter([0], [0], [0], marker="+", s=50)

        # http://stackoverflow.com/questions/11140163/python-matplotlib-plotting-a-3d-cube-a-sphere-and-a-vector
        # draw a vector
        from matplotlib.patches import FancyArrowPatch
        from mpl_toolkits.mplot3d import proj3d

        class Arrow3D(FancyArrowPatch):
            def __init__(self, xs, ys, zs, *args, **kwargs):
                FancyArrowPatch.__init__(self, (0, 0), (0, 0), *args, **kwargs)
                self._verts3d = xs, ys, zs

            def draw(self, renderer):
                xs3d, ys3d, zs3d = self._verts3d
                xs, ys, zs = proj3d.proj_transform(xs3d, ys3d, zs3d, renderer.M)
                self.set_positions((xs[0], ys[0]), (xs[1], ys[1]))
                FancyArrowPatch.draw(self, renderer)

        for v in self.candidate_basis_vectors:
            x, y, z = v.elems
            a = Arrow3D(
                [0, x],
                [0, y],
                [0, z],
                mutation_scale=10,
                lw=1,
                arrowstyle="-|>",
                color="k",
            )
            ax.add_artist(a)
            a = Arrow3D(
                [0, -x],
                [0, -y],
                [0, -z],
                mutation_scale=10,
                lw=1,
                arrowstyle="-|>",
                color="k",
            )
            ax.add_artist(a)

        x, y, z = zip(*self.candidate_basis_vectors)
        ax.scatter(x, y, z, marker=".", s=1)
        ax.scatter([-i for i in x], [-i for i in y], [-i for i in z], marker=".", s=1)
        pyplot.show()
