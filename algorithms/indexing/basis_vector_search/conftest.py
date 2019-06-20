from __future__ import absolute_import, division, print_function

import pytest
import random

from scitbx import matrix
from scitbx.math import euler_angles_as_matrix
from cctbx import sgtbx

random.seed(42)


def random_rotation(angle_min=0, angle_max=360):
    angles = [random.uniform(angle_min, angle_max) for i in range(3)]
    print("Rotation: ", angles)
    return euler_angles_as_matrix(angles, deg=True)


@pytest.fixture(params=["P2", "P3", "P6", "R3:h", "I23"])
def setup_rlp(request):
    random.seed(42)  # guaranteed to be random
    space_group = request.param
    # setup symmetry information
    sgi = sgtbx.space_group_info(symbol=space_group)
    cs = sgi.any_compatible_crystal_symmetry(volume=10000)
    cs = cs.best_cell()
    cs = cs.minimum_cell()
    uc = cs.unit_cell()

    B = matrix.sqr(uc.fractionalization_matrix()).transpose()
    U = random_rotation()
    A = U * B

    ms = cs.build_miller_set(d_min=2, anomalous_flag=True).expand_to_p1()
    rlp = A.elems * ms.indices().as_vec3_double()

    d = {}
    d["crystal_symmetry"] = cs
    d["rlp"] = rlp
    return d
