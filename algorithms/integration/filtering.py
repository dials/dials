#
# filtering.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
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
  ''' A class to do powder ring filtering. '''

  def __init__(self, unit_cell, space_group, d_min, width):
    ''' Initialise the filter. '''
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
    ''' True if within powder ring. '''
    from dials.array_family import flex
    result = flex.bool(len(d), False)
    for d_spacing in self.d_spacings:
      result = result | (flex.abs(d - d_spacing) < self.half_width)
    return result

  @classmethod
  def from_params(cls, params):
    ''' Factory method from phil. '''
    return PowderRingFilter(
      params.unit_cell,
      params.space_group.group(),
      params.d_min,
      params.width)


class MultiPowderRingFilter:
  ''' A class to encapsulate multiple powder ring filters '''

  def __init__(self):
    ''' Init the filter. '''
    self._filters = []

  def append(self, filter):
    ''' Add another powder ring filter. '''
    self._filters.append(filter)

  def __getitem__(self, index):
    ''' Get the powder ring filter at index. '''
    return self._filters[index]

  def __call__(self, d):
    ''' True if within powder ring. '''
    from dials.array_family import flex
    result = flex.bool(len(d), False)
    for filter in self:
      result = result | filter(d)
    return result

  def __len__(self):
    ''' The number of filters. '''
    return len(self._filters)

  def __iter__(self):
    ''' Iterate through filters. '''
    for i in range(len(self)):
      yield self[i]

  @classmethod
  def from_params(cls, params):
    ''' Factory method from phil. '''
    filters = cls()
    for i in range(len(params.powder.apply)):
      if params.powder.apply[i] == 'water_ice':
        filters.append(PowderRingFilter.from_params(params.powder.water_ice))
    return filters
