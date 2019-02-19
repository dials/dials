from __future__ import absolute_import, division

import pytest
import random

from scitbx import matrix
from scitbx.math import euler_angles_as_matrix
from cctbx import sgtbx

from dials.algorithms.indexing.basis_vector_search import strategies


def random_rotation(angle_min=0, angle_max=360):
    angles = [random.uniform(angle_min, angle_max) for i in xrange(3)]
    print("Rotation: ", angles)
    return euler_angles_as_matrix(angles, deg=True)


@pytest.fixture(params=["P2", "P3", "P6", "R3:h", "I23"])
def setup(request):
    space_group = request.param
    # setup symmetry information
    sgi = sgtbx.space_group_info(symbol=space_group)
    sg = sgi.group()
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


class TestStrategies(object):
    def check_results(self, unit_cell, basis_vectors):
        # check we have found a basis vector corresponding to each unit cell parameter
        for p in unit_cell.parameters()[:3]:
            found = False
            for v in basis_vectors:
                if p == pytest.approx(v.length(), abs=5e-1):
                    found = True
            assert found

    def test_fft1d(self, setup):
        max_cell = 1.3 * max(setup["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.fft1d(max_cell)
        basis_vectors = strategy.find_basis_vectors(setup["rlp"])
        self.check_results(setup["crystal_symmetry"].unit_cell(), basis_vectors)

    def test_fft3d(self, setup):
        max_cell = 1.3 * max(setup["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.fft3d(max_cell, n_points=256)
        basis_vectors = strategy.find_basis_vectors(setup["rlp"])
        self.check_results(setup["crystal_symmetry"].unit_cell(), basis_vectors)

    def test_real_space_grid_search(self, setup):
        max_cell = 1.3 * max(setup["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.real_space_grid_search(
            max_cell, target_unit_cell=setup["crystal_symmetry"].unit_cell()
        )
        basis_vectors = strategy.find_basis_vectors(setup["rlp"])
        self.check_results(setup["crystal_symmetry"].unit_cell(), basis_vectors)
