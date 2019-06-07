from __future__ import absolute_import, division

import pytest

from dials.algorithms.indexing.basis_vector_search import strategies


class TestStrategies(object):
    def check_results(self, unit_cell, basis_vectors):
        # check we have found a basis vector corresponding to each unit cell parameter
        for p in unit_cell.parameters()[:3]:
            found = False
            for v in basis_vectors:
                if p == pytest.approx(v.length(), abs=5e-1):
                    found = True
            assert found

    def test_fft1d(self, setup_rlp):
        max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.FFT1D(max_cell)
        basis_vectors, used = strategy.find_basis_vectors(setup_rlp["rlp"])
        self.check_results(setup_rlp["crystal_symmetry"].unit_cell(), basis_vectors)

    def test_fft3d(self, setup_rlp):
        max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.FFT3D(max_cell)
        basis_vectors, used = strategy.find_basis_vectors(setup_rlp["rlp"])
        self.check_results(setup_rlp["crystal_symmetry"].unit_cell(), basis_vectors)

    def test_real_space_grid_search(self, setup_rlp):
        max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
        strategy = strategies.RealSpaceGridSearch(
            max_cell, target_unit_cell=setup_rlp["crystal_symmetry"].unit_cell()
        )
        basis_vectors, used = strategy.find_basis_vectors(setup_rlp["rlp"])
        self.check_results(setup_rlp["crystal_symmetry"].unit_cell(), basis_vectors)
