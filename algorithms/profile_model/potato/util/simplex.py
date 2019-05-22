from __future__ import division

from scitbx import simplex


class SimpleSimplex(object):
    """
    Class to wrap some simplex stuff

    """

    def __init__(self, values, offset, evaluator, max_iter, tolerance=1e-10):
        """
        Init the simplex

        """
        self.n = len(values)
        self.x = values
        self.starting_simplex = self.generate_start(values, offset)
        optimizer = simplex.simplex_opt(
            dimension=self.n,
            matrix=self.starting_simplex,
            evaluator=evaluator,
            tolerance=tolerance,
            max_iter=max_iter,
        )
        self.x = optimizer.get_solution()

    def generate_start(self, values, offset):
        """
        Generate starting values

        """
        assert len(values) == len(offset)
        start = [values]
        for j, o in enumerate(offset):
            next = values.deep_copy()
            next[j] += o
            start.append(next)
        return start

    def get_solution(self):
        """
        Get the solution

        """
        return self.x
