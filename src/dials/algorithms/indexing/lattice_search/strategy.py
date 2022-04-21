"""Lattice search strategies."""


from __future__ import annotations


class Strategy:
    """A base class for lattice search strategies."""

    phil_help = None

    phil_scope = None

    def __init__(self, params=None, *args, **kwargs):
        """Construct the strategy.

        Args:
            params: an extracted PHIL scope containing the parameters
        """
        self._params = params
        if self._params is None and self.phil_scope is not None:
            self._params = self.phil_scope.extract()

    def find_crystal_models(self, reflections, experiments):
        """Find a list of likely crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

            experiments (dxtbx.model.experiment_list.ExperimentList):
                The experimental geometry models

        Returns:
            A list of candidate crystal models.
        """
        raise NotImplementedError()
