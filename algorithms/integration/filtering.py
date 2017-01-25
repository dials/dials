#
# filtering.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import absolute_import, division
from iotbx.phil import parse

# The phil scope
phil_scope = parse('''
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
      width = 0.06
        .type = float(value_min=0.0)
        .help = "The width of an ice ring (in d-spacing)."
    }

    apply = *none water_ice
      .type = choice(multi=True)
      .help = "The power ring filters to apply"
  }
''')


class PowderRingFilter:
  '''
  A class to do powder ring filtering.

  '''

  def __init__(self, unit_cell, space_group, d_min, width):
    '''
    Initialise the filter.

    :param unit_cell: The unit_cell of the powder rings
    :param space_group: The space group of the powder rings
    :param d_min: The maximum resolution to filter to
    :param width: The resolution width to filter around

    '''
    from cctbx.miller import index_generator
    from dials.array_family import flex
    assert(d_min > 0)
    assert(width > 0)

    # Correct unit cell
    unit_cell = space_group.average_unit_cell(unit_cell)

    # Generate a load of indices
    generator = index_generator(unit_cell, space_group.type(), False, d_min)
    indices = generator.to_array()

    # Compute d spacings and sort by resolution
    self.d_spacings = flex.sorted(unit_cell.d(indices))
    self.half_width = width / 2.0

  def __call__(self, d):
    '''
    True if within powder ring.

    :param d: The resolution
    :return: True/False in powder ring

    '''
    from dials.array_family import flex
    result = flex.bool(len(d), False)
    for d_spacing in self.d_spacings:
      result = result | (flex.abs(d - d_spacing) < self.half_width)
    return result

  @classmethod
  def from_params(cls, params):
    '''
    Factory method from phil.

    :param params: The input phil parameters
    :return: The powder ring filter

    '''
    return PowderRingFilter(
      params.unit_cell,
      params.space_group.group(),
      params.d_min,
      params.width)


class MultiPowderRingFilter:
  '''
  A class to encapsulate multiple powder ring filters

  '''

  def __init__(self):
    '''
    Init the filter.

    '''
    self._filters = []

  def append(self, filter):
    '''
    Add another powder ring filter.

    :param filter: The filter to add

    '''
    self._filters.append(filter)

  def __getitem__(self, index):
    '''
    Get the powder ring filter at index.

    :param index: The index of the filter
    :return: The requested filter

    '''
    return self._filters[index]

  def __call__(self, d):
    '''
    True if within powder ring.

    :param d: The resolution
    :return: True/False if within a powder ring

    '''
    from dials.array_family import flex
    result = flex.bool(len(d), False)
    for filter in self:
      result = result | filter(d)
    return result

  def __len__(self):
    '''
    :return: The number of filters.

    '''
    return len(self._filters)

  def __iter__(self):
    '''
    Iterate through filters.

    '''
    for i in range(len(self)):
      yield self[i]

  @classmethod
  def from_params(cls, params):
    '''
    Factory method from phil.

    :param params: The input phil parameters
    :return: The powder ring filter

    '''
    filters = cls()
    for i in range(len(params.powder.apply)):
      if params.powder.apply[i] == 'water_ice':
        filters.append(PowderRingFilter.from_params(params.powder.water_ice))
    return filters


class IceRingFilter:
  '''
  A class to do ice ring filtering

  '''

  def __init__(self):
    '''
    Initialise the filter.

    :param width: The resolution width to filter around

    '''
    # Hexagonal ice ring resolution ranges in 1/d^2
    self.ice_rings = [
      (0.0640, 0.0690),
      (0.0710, 0.0780),
      (0.0825, 0.0880),
      ( 0.138,  0.144),
      ( 0.190,  0.205),
      ( 0.228,  0.240),
      ( 0.262,  0.266),
      ( 0.267,  0.278),
      ( 0.280,  0.288),
      ( 0.337,  0.341),
      ( 0.429,  0.435),
      ( 0.459,  0.466),
      ( 0.478,  0.486),
      ( 0.531,  0.537),
    ]

  def __call__(self, d):
    '''
    True if within powder ring.

    :param d: The resolution
    :return: True/False in powder ring

    '''
    from dials.array_family import flex
    result = flex.bool(len(d), False)
    d2 = 1.0 / d**2
    for ice_ring in self.ice_rings:
      result = result | (d2 >= ice_ring[0]) & (d2 <= ice_ring[1])
    return result
