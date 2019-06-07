from __future__ import absolute_import, division, print_function

import warnings

import dials_algorithms_indexing_ext as ext

sampling_volume_map = ext.sampling_volume_map


class DialsIndexError(RuntimeError):
    pass


class DialsIndexRefineError(DialsIndexError):
    pass


def index_reflections(reflections, experiments, d_min=None, tolerance=0.3):
    from dials.algorithms.indexing import assign_indices

    warnings.warn(
        "index_reflections is deprectated, use "
        "dials.algorithms.indexing.assign_indices.AssignIndicesGlobal instead",
        DeprecationWarning,
        stacklevel=2,
    )
    index = assign_indices.AssignIndicesGlobal(tolerance=tolerance)
    index(reflections, experiments)
