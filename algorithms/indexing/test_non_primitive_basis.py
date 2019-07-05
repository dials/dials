from __future__ import absolute_import, division, print_function

import pytest

from cctbx import sgtbx
from cctbx.sgtbx import bravais_types
import scitbx.matrix
from dxtbx.model import Crystal, Experiment, ExperimentList
from dials.algorithms.indexing import non_primitive_basis
from dials.algorithms.indexing import assign_indices
from dials.array_family import flex


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_detect(space_group_symbol):
    sgi = sgtbx.space_group_info(space_group_symbol)
    cs = sgi.any_compatible_crystal_symmetry(volume=1000)
    ms = cs.build_miller_set(anomalous_flag=True, d_min=1).expand_to_p1()
    result = non_primitive_basis.detect(ms.indices())
    if sgi.group().conventional_centring_type_symbol() != "P":
        assert result is not None
        assert isinstance(result, scitbx.matrix.sqr)
        assert result.n == (3, 3)
    else:
        assert result is None


@pytest.mark.parametrize("space_group_symbol", bravais_types.acentric)
def test_correct(space_group_symbol):

    sgi = sgtbx.space_group_info(space_group_symbol)
    cs = sgi.any_compatible_crystal_symmetry(volume=1000)
    ms = cs.build_miller_set(anomalous_flag=True, d_min=1).expand_to_p1()

    # the reciprocal matrix
    B = scitbx.matrix.sqr(cs.unit_cell().fractionalization_matrix()).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    expts = ExperimentList([Experiment(crystal=crystal)])

    refl = flex.reflection_table()
    refl["miller_index"] = ms.indices()
    refl["rlp"] = B.elems * ms.indices().as_vec3_double()
    refl["imageset_id"] = flex.int(len(refl))
    refl["xyzobs.mm.value"] = flex.vec3_double(len(refl))

    non_primitive_basis.correct(expts, refl, assign_indices.AssignIndicesGlobal())

    cs_corrected = expts.crystals()[0].get_crystal_symmetry()
    assert cs_corrected.change_of_basis_op_to_primitive_setting().is_identity_op()
    assert (
        cs.change_of_basis_op_to_primitive_setting().apply(ms.indices())
        == refl["miller_index"]
    )
