"""
Tools for analysis of missing reflections.
"""

from __future__ import annotations

from annlib_ext import AnnAdaptor

import cctbx.miller
from boost_adaptbx import graph
from boost_adaptbx.graph import connected_component_algorithm as cca

from dials.array_family import flex


def connected_components(
    miller_array: cctbx.miller.array,
) -> [{}]:
    """
    Identify connected regions of missing reflections in the asymmetric unit.

    This is achieved by first generating the complete set of possible miller indices,
    then performing connected components analysis on a graph of nearest neighbours in
    the list of missing reflections.

    Args:
        miller_array:  The input list of reflections.

    Returns:
        The list of miller sets for each connected region of missing reflections. The
        first item in the list will be the complete set of all possible miller indices.
    """

    # Map to primitive setting for centred cells, otherwise true missing reflections
    # won't be identified as connected as a result of being separated by systematically
    # absent reflections.
    cb_op_to_primitive = miller_array.change_of_basis_op_to_primitive_setting()
    miller_array = miller_array.change_basis(cb_op_to_primitive)

    # First generate the missing_set of reflections. We want the full sphere of missing
    # reflections to allow us to find connected regions that cross the boundary of the
    # asu.
    unique = miller_array.unique_under_symmetry().map_to_asu()
    unique = unique.generate_bijvoet_mates()
    complete_set = unique.complete_set()
    missing_set = complete_set.lone_set(unique)
    missing_set = missing_set.expand_to_p1().customized_copy(
        crystal_symmetry=missing_set.crystal_symmetry()
    )

    if missing_set.size() == 0:
        return complete_set, []

    # Now find the nearest neighbours.
    mi = missing_set.indices().as_vec3_double().as_double()
    k = 6
    ann = AnnAdaptor(data=mi, dim=3, k=k)
    ann.query(mi)

    # Construct the graph of connected missing reflections
    g = graph.adjacency_list(
        graph_type="undirected",
        vertex_type="vector",
        edge_type="set",
    )
    distance_cutoff = 2**0.5
    for i in range(missing_set.size()):
        ik = i * k
        for i_ann in range(k):
            if ann.distances[ik + i_ann] <= distance_cutoff:
                j = ann.nn[ik + i_ann]
                g.add_edge(i, j)

    # Now do the connected components analysis, filtering out lone missing reflections
    components = [c for c in cca.connected_components(graph=g) if len(c) > 1]

    # Determine the unique miller indices for each component within the asu
    unique_mi = []
    unique_ms = []
    for i, c in enumerate(components):
        ms = (
            missing_set.select(flex.size_t(list(c)))
            .customized_copy(crystal_symmetry=miller_array)
            .as_non_anomalous_set()
            .map_to_asu()
        )
        ms = ms.unique_under_symmetry()
        mi = set(ms.indices())
        if mi not in unique_mi:
            unique_ms.append(ms)
            unique_mi.append(mi)

    # Sort connected regions by size
    unique_ms = sorted(unique_ms, key=lambda ms: ms.size(), reverse=True)

    # Map indices back to input setting
    cb_op_primitive_inp = cb_op_to_primitive.inverse()
    return (
        unique.as_non_anomalous_set().complete_set().change_basis(cb_op_primitive_inp),
        [ms.change_basis(cb_op_primitive_inp) for ms in unique_ms],
    )
