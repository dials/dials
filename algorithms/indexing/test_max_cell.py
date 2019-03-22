from __future__ import absolute_import, division, print_function

import pytest
import random

from cctbx import sgtbx
from cctbx.sgtbx import bravais_types
import scitbx.matrix

from dxtbx.model import Crystal, Experiment, ExperimentList
from dials.array_family import flex
from dials.algorithms.indexing.max_cell import find_max_cell


@pytest.fixture(params=bravais_types.acentric)
def setup(request):
    space_group_symbol = request.param
    sgi = sgtbx.space_group_info(space_group_symbol)
    cs = sgi.any_compatible_crystal_symmetry(volume=random.randint(1e4, 1e6))
    ms = cs.build_miller_set(anomalous_flag=True, d_min=3).expand_to_p1()

    # the reciprocal matrix
    B = scitbx.matrix.sqr(cs.unit_cell().fractionalization_matrix()).transpose()
    crystal = Crystal(B, sgtbx.space_group())
    expts = ExperimentList([Experiment(crystal=crystal)])

    # randomly select 25% of reflections
    ms = ms.select(flex.random_permutation(ms.size())[: int(0.25 * ms.size())])

    refl = flex.reflection_table()
    refl["rlp"] = B.elems * ms.indices().as_vec3_double()
    refl["imageset_id"] = flex.int(len(refl))
    refl["xyzobs.mm.value"] = flex.vec3_double(len(refl))

    d = {}
    d["crystal_symmetry"] = cs
    d["reflections"] = refl
    return d


def test_max_cell(setup):
    reflections = setup["reflections"]
    crystal_symmetry = setup["crystal_symmetry"]

    max_cell_multiplier = 1.3
    max_cell = find_max_cell(reflections, max_cell_multiplier=max_cell_multiplier)

    known_max_cell = max(
        crystal_symmetry.primitive_setting().unit_cell().parameters()[:3]
    )
    assert max_cell.max_cell > known_max_cell
