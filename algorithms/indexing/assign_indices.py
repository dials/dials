from __future__ import absolute_import, division

import abc

from cctbx.array_family import flex
import dials_algorithms_indexing_ext as ext

from dials.algorithms.indexing import DialsIndexError


class AssignIndicesStrategy(object):

    __metaclass__ = abc.ABCMeta

    def __init__(self, d_min=None):
        self._d_min = d_min

    @abc.abstractmethod
    def __call__(self, reciprocal_lattice_vectors):
        pass


class AssignIndicesGlobal(AssignIndicesStrategy):
    def __init__(self, tolerance=0.3):
        super(AssignIndicesGlobal, self).__init__()
        self._tolerance = tolerance

    def __call__(self, reflections, experiments, d_min=None):
        reciprocal_lattice_points = reflections["rlp"]
        reflections["miller_index"] = flex.miller_index(len(reflections), (0, 0, 0))
        if d_min is not None:
            d_spacings = 1 / reciprocal_lattice_points.norms()
            inside_resolution_limit = d_spacings > d_min
        else:
            inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
        sel = inside_resolution_limit & (reflections["id"] == -1)
        isel = sel.iselection()
        rlps = reciprocal_lattice_points.select(isel)
        refs = reflections.select(isel)
        phi = refs["xyzobs.mm.value"].parts()[2]

        UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])
        imgset_ids = reflections["imageset_id"].select(sel)

        for i_imgset, imgset in enumerate(experiments.imagesets()):
            sel_imgset = imgset_ids == i_imgset

            result = ext.AssignIndices(
                rlps.select(sel_imgset),
                phi.select(sel_imgset),
                UB_matrices,
                tolerance=self._tolerance,
            )

            miller_indices = result.miller_indices()
            crystal_ids = result.crystal_ids()

            expt_ids = flex.int(crystal_ids.size(), -1)
            for i_cryst, cryst in enumerate(experiments.crystals()):
                sel_cryst = crystal_ids == i_cryst
                for i_expt in experiments.where(crystal=cryst, imageset=imgset):
                    expt_ids.set_selected(sel_cryst, i_expt)

            reflections["miller_index"].set_selected(
                isel.select(sel_imgset), miller_indices
            )
            reflections["id"].set_selected(isel.select(sel_imgset), expt_ids)
            reflections.set_flags(
                reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
            )
            reflections["id"].set_selected(reflections["miller_index"] == (0, 0, 0), -1)


class AssignIndicesLocal(AssignIndicesStrategy):
    def __init__(
        self, d_min=None, epsilon=0.05, delta=8, l_min=0.8, nearest_neighbours=20
    ):
        super(AssignIndicesLocal, self).__init__()
        self._epsilon = epsilon
        self._delta = delta
        self._l_min = l_min
        self._nearest_neighbours = nearest_neighbours

    def __call__(self, reflections, experiments, d_min=None):
        from scitbx import matrix
        from libtbx.math_utils import nearest_integer as nint

        reciprocal_lattice_points = reflections["rlp"]
        if "miller_index" not in reflections:
            reflections["miller_index"] = flex.miller_index(len(reflections))
        if d_min is not None:
            d_spacings = 1 / reciprocal_lattice_points.norms()
            inside_resolution_limit = d_spacings > d_min
        else:
            inside_resolution_limit = flex.bool(reciprocal_lattice_points.size(), True)
        sel = inside_resolution_limit & (reflections["id"] == -1)
        isel = sel.iselection()
        rlps = reciprocal_lattice_points.select(isel)
        refs = reflections.select(isel)
        phi = refs["xyzobs.mm.value"].parts()[2]

        if len(rlps) <= self._nearest_neighbours:
            raise DialsIndexError(
                "index_assignment.local.nearest_neighbour must be smaller than the number of accepted reflections (%d)"
                % len(rlps)
            )

        UB_matrices = flex.mat3_double([cm.get_A() for cm in experiments.crystals()])

        result = ext.AssignIndicesLocal(
            rlps,
            phi,
            UB_matrices,
            epsilon=self._epsilon,
            delta=self._delta,
            l_min=self._l_min,
            nearest_neighbours=self._nearest_neighbours,
        )
        miller_indices = result.miller_indices()
        crystal_ids = result.crystal_ids()
        hkl = miller_indices.as_vec3_double().iround()

        assert miller_indices.select(crystal_ids < 0).all_eq((0, 0, 0))

        for i_cryst in set(crystal_ids):
            if i_cryst < 0:
                continue

            A = matrix.sqr(experiments[i_cryst].crystal.get_A())
            A_inv = A.inverse()

            cryst_sel = crystal_ids == i_cryst
            rlp_sel = rlps.select(cryst_sel)
            hkl_sel = hkl.select(cryst_sel).as_vec3_double()

            d_sel = 1 / rlp_sel.norms()
            d_perm = flex.sort_permutation(d_sel, reverse=True)

            hf_0 = A_inv * rlp_sel[d_perm[0]]
            h_0 = matrix.col([nint(j) for j in hf_0.elems])
            offset = h_0 - matrix.col(hkl_sel[d_perm[0]])
            # print "offset:", offset.elems

            h = hkl_sel + flex.vec3_double(hkl_sel.size(), offset.elems)

            refs["miller_index"].set_selected(
                cryst_sel, flex.miller_index(list(h.iround()))
            )
            refs["id"].set_selected(cryst_sel, i_cryst)

        crystal_ids.set_selected(crystal_ids < 0, -1)
        refs["id"] = crystal_ids
        refs["miller_index"].set_selected(crystal_ids < 0, (0, 0, 0))

        reflections["miller_index"].set_selected(isel, refs["miller_index"])
        reflections["id"].set_selected(isel, refs["id"])
        reflections.set_flags(
            reflections["miller_index"] != (0, 0, 0), reflections.flags.indexed
        )
