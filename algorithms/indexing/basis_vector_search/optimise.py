from __future__ import absolute_import, division
from __future__ import print_function

import math

from scitbx.array_family import flex
from scitbx import lbfgs


def optimise_basis_vectors(reciprocal_lattice_points, vectors):
    optimised = flex.vec3_double()
    for vector in vectors:
        minimised = BasisVectorMinimiser(reciprocal_lattice_points, vector)
        optimised.append(tuple(minimised.x))
    functionals = flex.double(minimised.target.compute_functional(v) for v in vectors)
    perm = flex.sort_permutation(functionals)
    optimised = optimised.select(perm)
    return optimised


# Optimise the initial basis vectors as per equation 11.4.3.4 of
# Otwinowski et al, International Tables Vol. F, chapter 11.4 pp. 282-295
class BasisVectorTarget(object):
    def __init__(self, reciprocal_lattice_points):
        self.reciprocal_lattice_points = reciprocal_lattice_points
        self._xyz_parts = self.reciprocal_lattice_points.parts()

    def compute_functional(self, vector):
        two_pi_S_dot_v = 2 * math.pi * self.reciprocal_lattice_points.dot(vector)
        return -flex.sum(flex.cos(two_pi_S_dot_v))

    def compute_functional_and_gradients(self, vector):
        assert len(vector) == 3
        two_pi_S_dot_v = 2 * math.pi * self.reciprocal_lattice_points.dot(vector)
        f = -flex.sum(flex.cos(two_pi_S_dot_v))
        sin_part = flex.sin(two_pi_S_dot_v)
        g = flex.double(
            [flex.sum(2 * math.pi * self._xyz_parts[i] * sin_part) for i in range(3)]
        )
        return f, g


class BasisVectorMinimiser(object):
    def __init__(
        self,
        reciprocal_lattice_points,
        vector,
        lbfgs_termination_params=None,
        lbfgs_core_params=lbfgs.core_parameters(m=20),
    ):
        self.reciprocal_lattice_points = reciprocal_lattice_points
        if not isinstance(vector, flex.double):
            self.x = flex.double(vector)
        else:
            self.x = vector.deep_copy()
        self.n = len(self.x)
        assert self.n == 3
        self.target = BasisVectorTarget(self.reciprocal_lattice_points)
        self.minimizer = lbfgs.run(
            target_evaluator=self,
            termination_params=lbfgs_termination_params,
            core_params=lbfgs_core_params,
        )

    def compute_functional_and_gradients(self):
        f, g = self.target.compute_functional_and_gradients(tuple(self.x))
        return f, g
