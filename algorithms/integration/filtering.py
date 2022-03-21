from __future__ import annotations

from cctbx import uctbx
from cctbx.miller import index_generator
from iotbx.phil import parse

from dials.array_family import flex

# The phil scope
phil_scope = parse(
    """
  powder {

    water_ice {
      unit_cell = 4.498,4.498,7.338,90,90,120
        .type = unit_cell
        .help = "The unit cell to generate d_spacings for ice rings."
      space_group = 194
        .type = space_group
        .help = "The space group used to generate d_spacings for ice rings."
      d_min = 1
        .type = float(value_min=0.0)
        .help = "The minimum resolution to filter ice rings"
      width = 0.002
        .type = float(value_min=0.0)
        .help = "The width of an ice ring (in 1/d^2)."
    }

    apply = *none water_ice
      .type = choice(multi=True)
      .help = "The power ring filters to apply"
  }
"""
)


class PowderRingFilter:
    """
    A class to do powder ring filtering.
    """

    def __init__(self, unit_cell, space_group, d_min, width):
        """
        Initialise the filter.

        :param unit_cell: The unit_cell of the powder rings
        :param space_group: The space group of the powder rings
        :param d_min: The maximum resolution to filter to
        :param width: The resolution width to filter around
        """
        assert d_min > 0
        assert width > 0

        # Correct unit cell
        unit_cell = space_group.average_unit_cell(unit_cell)

        self.half_width = width / 2.0
        d_min = uctbx.d_star_sq_as_d(uctbx.d_as_d_star_sq(d_min) + self.half_width)

        # Generate a load of indices
        generator = index_generator(unit_cell, space_group.type(), False, d_min)
        indices = generator.to_array()

        # Compute d spacings and sort by resolution
        self.d_star_sq = flex.sorted(unit_cell.d_star_sq(indices))

    def __call__(self, d):
        """
        True if within powder ring.

        :param d: The resolution
        :return: True/False in powder ring
        """
        result = flex.bool(len(d), False)
        d_star_sq = uctbx.d_as_d_star_sq(d)
        for ds2 in self.d_star_sq:
            result = result | (flex.abs(d_star_sq - ds2) < self.half_width)
        return result


class IceRingFilter:
    """
    A class to do ice ring filtering
    """

    def __init__(self):
        """
        Initialise the filter.

        :param width: The resolution width to filter around
        """
        # Hexagonal ice ring resolution ranges in 1/d^2
        self.ice_rings = [
            (0.0640, 0.0690),
            (0.0710, 0.0780),
            (0.0825, 0.0880),
            (0.138, 0.144),
            (0.190, 0.205),
            (0.228, 0.240),
            (0.262, 0.266),
            (0.267, 0.278),
            (0.280, 0.288),
            (0.337, 0.341),
            (0.429, 0.435),
            (0.459, 0.466),
            (0.478, 0.486),
            (0.531, 0.537),
        ]

    def __call__(self, d):
        """
        True if within powder ring.

        :param d: The resolution
        :return: True/False in powder ring
        """
        result = flex.bool(len(d), False)
        d2 = 1.0 / flex.pow2(d)
        for ice_ring in self.ice_rings:
            result = result | (d2 >= ice_ring[0]) & (d2 <= ice_ring[1])
        return result
