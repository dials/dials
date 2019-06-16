from __future__ import absolute_import, division, print_function

import pytest
import random

from dials.algorithms.indexing.basis_vector_search import optimise, strategies


def test_optimise_basis_vectors(setup_rlp):
    random.seed(42)  # guaranteed to be random
    max_cell = 1.3 * max(setup_rlp["crystal_symmetry"].unit_cell().parameters()[:3])
    rlp = setup_rlp["rlp"]
    strategy = strategies.FFT3D(max_cell, n_points=256)
    basis_vectors, used = strategy.find_basis_vectors(rlp)

    for v in basis_vectors:
        target = optimise.BasisVectorTarget(rlp)
        f, g = target.compute_functional_and_gradients(v)
        assert f == target.compute_functional(v)
        g_fd = _gradient_fd(target, v)
        assert list(g) == pytest.approx(g_fd, rel=1)

        minimised = optimise.BasisVectorMinimiser(rlp, v)
        # check that the minimised vectors are broadly similar to the starting vectors
        assert v.length() == pytest.approx(minimised.x.norm(), rel=2e-2)
        assert v.elems == pytest.approx(minimised.x, abs=1.2)

    optimised = optimise.optimise_basis_vectors(rlp, basis_vectors)
    assert len(optimised) == len(basis_vectors)


def _gradient_fd(target, vector, eps=1e-6):
    grads = []
    for i in range(len(vector)):
        v = list(vector)
        v[i] -= eps
        tm, _ = target.compute_functional_and_gradients(v)
        v[i] += 2 * eps
        tp, _ = target.compute_functional_and_gradients(v)
        grads.append((tp - tm) / (2 * eps))
    return grads
