from __future__ import absolute_import, division, print_function

import abc
import copy
import logging
import math

import libtbx
from libtbx import group_args
from scitbx import matrix
from scitbx.array_family import flex
from dxtbx.model import Crystal
from dials.algorithms.indexing.compare_orientation_matrices import (
    difference_rotation_matrix_axis_angle,
)

logger = logging.getLogger(__name__)


class Result(group_args):
    pass


def filter_doubled_cell(solutions):
    accepted_solutions = []
    for i1, s1 in enumerate(solutions):
        doubled_cell = False
        for (m1, m2, m3) in (
            (2, 1, 1),
            (1, 2, 1),
            (1, 1, 2),
            (2, 2, 1),
            (2, 1, 2),
            (1, 2, 2),
            (2, 2, 2),
        ):
            if doubled_cell:
                break
            a, b, c = (matrix.col(v) for v in s1.crystal.get_real_space_vectors())
            new_cryst = Crystal(
                real_space_a=1 / m1 * a,
                real_space_b=1 / m2 * b,
                real_space_c=1 / m3 * c,
                space_group=s1.crystal.get_space_group(),
            )
            new_unit_cell = new_cryst.get_unit_cell()
            for s2 in solutions:
                if s2 is s1:
                    continue
                if new_unit_cell.is_similar_to(
                    s2.crystal.get_unit_cell(), relative_length_tolerance=0.05
                ):
                    R, axis, angle, cb = difference_rotation_matrix_axis_angle(
                        new_cryst, s2.crystal
                    )
                    if (angle < 1) and (s1.n_indexed < (1.1 * s2.n_indexed)):
                        doubled_cell = True
                        break

        if not doubled_cell:
            accepted_solutions.append(s1)

    return accepted_solutions


class ModelRank(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self):
        self.all_solutions = []

    def append(self, item):
        self.all_solutions.append(item)

    def extend(self, items):
        self.all_solutions.extend(items)

    @abc.abstractmethod
    def best_model(self):
        pass

    @abc.abstractmethod
    def __str__(self):
        pass


# Tracker for solutions based on code in rstbx/dps_core/basis_choice.py
class ModelRankFilter(ModelRank):
    def __init__(
        self,
        check_doubled_cell=True,
        likelihood_cutoff=0.8,
        volume_cutoff=1.25,
        n_indexed_cutoff=0.9,
    ):
        super(ModelRankFilter, self).__init__()
        self.check_doubled_cell = check_doubled_cell
        self.likelihood_cutoff = likelihood_cutoff
        self.volume_cutoff = volume_cutoff
        self.n_indexed_cutoff = n_indexed_cutoff
        self.filtered_solutions = []

    def append(self, item):
        super(ModelRankFilter, self).append(item)
        self.update_analysis()

    def extend(self, items):
        super(ModelRankFilter, self).extend(items)
        self.update_analysis()

    def __len__(self):
        return len(self.filtered_solutions)

    def filter_by_likelihood(self, solutions):
        best_likelihood = max(s.model_likelihood for s in solutions)
        offset = 0
        while (best_likelihood + offset) <= 0:
            offset += 1
        return [
            s
            for s in solutions
            if (s.model_likelihood + offset)
            >= (self.likelihood_cutoff * (best_likelihood + offset))
        ]

    def filter_by_volume(self, solutions):
        # filter by volume - prefer solutions with a smaller unit cell
        min_volume = min(s.crystal.get_unit_cell().volume() for s in solutions)
        return [
            s
            for s in solutions
            if s.crystal.get_unit_cell().volume() < (self.volume_cutoff * min_volume)
        ]

    def filter_by_n_indexed(self, solutions, n_indexed_cutoff=None):
        if n_indexed_cutoff is None:
            n_indexed_cutoff = self.n_indexed_cutoff
        # filter by number of indexed reflections - prefer solutions that
        # account for more of the diffracted spots
        max_n_indexed = max(s.n_indexed for s in solutions)
        return [s for s in solutions if s.n_indexed >= n_indexed_cutoff * max_n_indexed]

    def update_analysis(self):
        # pre-filter out solutions that only account for a very small
        # percentage of the indexed spots relative to the best one
        self.filtered_solutions = self.filter_by_n_indexed(
            self.all_solutions, n_indexed_cutoff=0.05
        )  # 5 percent

        if self.check_doubled_cell:
            self.filtered_solutions = filter_doubled_cell(self.filtered_solutions)

        self.filtered_solutions = self.filter_by_likelihood(self.filtered_solutions)

        self.filtered_solutions = self.filter_by_volume(self.filtered_solutions)

        self.filtered_solutions = self.filter_by_n_indexed(self.filtered_solutions)

        return

    def best_model(self):
        self.best_filtered_liklihood = max(
            s.model_likelihood for s in self.filtered_solutions
        )

        solutions = [
            s
            for s in self.filtered_solutions
            if s.model_likelihood == self.best_filtered_liklihood
        ]
        return solutions[0]

    def __str__(self):
        rows = []
        rows.append(
            ["unit_cell", "volume", "n_indexed", "fraction_indexed", "likelihood"]
        )

        for i, s in enumerate(self.all_solutions):
            s = self.all_solutions[i]
            rows.append(
                [
                    format(
                        s.crystal.get_unit_cell(),
                        "{:.2f} {:.2f} {:.2f} {:.1f} {:.1f} {:.1f}",
                    ),
                    "%.0f" % s.crystal.get_unit_cell().volume(),
                    str(s.n_indexed),
                    "%.0f" % (s.fraction_indexed * 100),
                    "%.2f" % s.model_likelihood,
                ]
            )

        from libtbx import table_utils

        return table_utils.format(rows=rows, has_header=True)


class ModelRankWeighted(ModelRank):
    def __init__(self, power=2, volume_weight=1, n_indexed_weight=1, rmsd_weight=1):
        super(ModelRankWeighted, self).__init__()
        self.volume_weight = volume_weight
        self.n_indexed_weight = n_indexed_weight
        self.rmsd_weight = rmsd_weight
        self.power = power

    def __len__(self):
        return len(self.all_solutions)

    def score_by_volume(self, reverse=False):
        # smaller volume = better
        volumes = flex.double(
            s.crystal.get_unit_cell().volume() for s in self.all_solutions
        )
        score = flex.log(volumes) / math.log(2)
        return self.volume_weight * (score - flex.min(score))

    def score_by_rmsd_xy(self, reverse=False):
        # smaller rmsds = better
        rmsd_x, rmsd_y, rmsd_z = flex.vec3_double(
            s.rmsds for s in self.all_solutions
        ).parts()
        rmsd_xy = flex.sqrt(flex.pow2(rmsd_x) + flex.pow2(rmsd_y))
        score = flex.log(rmsd_xy) / math.log(2)
        return self.rmsd_weight * (score - flex.min(score))

    def score_by_fraction_indexed(self, reverse=False):
        # more indexed reflections = better
        fraction_indexed = flex.double(s.fraction_indexed for s in self.all_solutions)
        score = flex.log(fraction_indexed) / math.log(2)
        return self.n_indexed_weight * (-score + flex.max(score))

    def best_model(self):
        scores = self.combined_scores()
        perm = flex.sort_permutation(scores)
        return self.all_solutions[perm[0]]

    def combined_scores(self):
        scores = sum(
            flex.pow(score.as_double(), self.power)
            for score in (
                self.score_by_fraction_indexed(),
                self.score_by_volume(),
                self.score_by_rmsd_xy(),
            )
        )
        return scores

    def __str__(self):
        rows = []
        rows.append(
            [
                "unit_cell",
                "volume",
                "volume score",
                "#indexed",
                "% indexed",
                "% indexed score",
                "rmsd_xy",
                "rmsd_xy score",
                "overall score",
            ]
        )

        score_by_fraction_indexed = self.score_by_fraction_indexed()
        score_by_volume = self.score_by_volume()
        score_by_rmsd_xy = self.score_by_rmsd_xy()
        combined_scores = self.combined_scores()

        perm = flex.sort_permutation(combined_scores)

        rmsd_x, rmsd_y, rmsd_z = flex.vec3_double(
            s.rmsds for s in self.all_solutions
        ).parts()
        rmsd_xy = flex.sqrt(flex.pow2(rmsd_x) + flex.pow2(rmsd_y))

        for i in perm:
            s = self.all_solutions[i]
            rows.append(
                [
                    format(
                        s.crystal.get_unit_cell(),
                        "{:.2f} {:.2f} {:.2f} {:.1f} {:.1f} {:.1f}",
                    ),
                    "%.0f" % s.crystal.get_unit_cell().volume(),
                    "%.2f" % score_by_volume[i],
                    str(s.n_indexed),
                    "%.0f" % (s.fraction_indexed * 100),
                    "%.2f" % score_by_fraction_indexed[i],
                    "%.2f" % rmsd_xy[i],
                    "%.2f" % score_by_rmsd_xy[i],
                    "%.2f" % combined_scores[i],
                ]
            )

        from libtbx import table_utils

        return table_utils.format(rows=rows, has_header=True)


class Strategy(object):

    __metaclass__ = abc.ABCMeta

    @abc.abstractmethod
    def evaluate(self, experiments, reflections):
        pass


class ModelEvaluation(Strategy):
    def __init__(self, refinement_params):
        self._params = copy.deepcopy(refinement_params)
        # override several parameters, mainly for speed
        self._params.refinement.parameterisation.auto_reduction.action = "fix"
        self._params.refinement.parameterisation.scan_varying = False
        self._params.refinement.refinery.max_iterations = 4
        self._params.refinement.reflections.reflections_per_degree = min(
            self._params.refinement.reflections.reflections_per_degree, 20
        )
        if self._params.refinement.reflections.outlier.block_width is libtbx.Auto:
            # auto block_width determination is potentially too expensive to do at
            # this stage: instead set separate_blocks=False and increase value
            # of tukey.iqr_multiplier to be more tolerant of outliers
            self._params.refinement.reflections.outlier.separate_blocks = False
            self._params.refinement.reflections.outlier.tukey.iqr_multiplier = (
                2 * self._params.refinement.reflections.outlier.tukey.iqr_multiplier
            )

    def evaluate(self, experiments, reflections):

        from dials.algorithms.refinement import RefinerFactory

        indexed_reflections = reflections.select(reflections["id"] > -1)
        reflogger = logging.getLogger("dials.algorithms.refinement")
        level = reflogger.getEffectiveLevel()
        reflogger.setLevel(logging.ERROR)
        try:
            refiner = RefinerFactory.from_parameters_data_experiments(
                self._params, indexed_reflections, experiments
            )
            refiner.run()
        except (RuntimeError, ValueError):
            return
        else:
            rmsds = refiner.rmsds()
            xy_rmsds = math.sqrt(rmsds[0] ** 2 + rmsds[1] ** 2)
            model_likelihood = 1.0 - xy_rmsds
            result = Result(
                model_likelihood=model_likelihood,
                crystal=experiments.crystals()[0],
                rmsds=rmsds,
                n_indexed=len(indexed_reflections),
                fraction_indexed=float(len(indexed_reflections)) / len(reflections),
                hkl_offset=(0, 0, 0),
            )
            return result
        finally:
            reflogger.setLevel(level)
