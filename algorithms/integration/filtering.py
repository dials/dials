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
    generator = index_generator(unit_cell, space_group.type(), false, d_min)
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
    return filters_[index]

  def __call__(self, d):
    ''' True if within powder ring. '''
    result = flex.bool(len(d), False)
    for filter in self:
      result = result | filter(d)
    return result

  def __len__(self):
    ''' The number of filters. '''
    return len(self.filters)

  def __iter__(self):
    ''' Iterate through filters. '''
    for i in range(len(self)):
      yield self[i]
