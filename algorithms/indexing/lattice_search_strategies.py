"""Lattice search strategies."""

from __future__ import absolute_import, division
from __future__ import print_function

import abc
import logging

from libtbx import phil

logger = logging.getLogger(__name__)


class Strategy(object):
    """A base class for lattice search strategies."""

    __metaclass__ = abc.ABCMeta

    phil_scope = None

    def __init__(self, params=None, *args, **kwargs):
        """Construct the strategy.

        Args:
            params: an extracted PHIL scope containing the parameters

        """
        self._params = params
        if self._params is None and self.phil_scope is not None:
            self._params = self.phil_scope.extract()

    @abc.abstractmethod
    def find_crystal_models(self, reflections):
        """Find a list of likely crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

        Returns:
            A list of candidate crystal models.

        """
        pass


low_res_spot_match_phil_str = """\

"""


class LowResSpotMatch(Strategy):
    """Lattice search by matching low resolution spots to candidate indices
    based on resolution and reciprocal space distance between observed spots.
    """

    phil_scope = phil.parse(low_res_spot_match_phil_str)

    def __init__(self, target_symmetry_primitive, params=None, *args, **kwargs):
        """Construct a real_space_grid_search object.

        Args:
            target_symmetry_primitive (cctbx.crystal.symmetry): The target
                crystal symmetry and unit cell

        """
        super(LowResSpotMatch, self).__init__(params=params, *args, **kwargs)
        self._target_symmetry_primitive = target_symmetry_primitive

    def find_crystal_models(self, reflections):
        """Find a list of candidate crystal models.

        Args:
            reflections (dials.array_family.flex.reflection_table):
                The found spots centroids and associated data

        """

        return []
