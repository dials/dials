"""
Functions to help with calculating batch properties for experiments objects.
"""
from math import ceil, floor, log
import logging
import copy as copy
from dials.array_family import flex
from dials.util import Sorry

logger = logging.getLogger('dials')

def assign_batches_to_reflections(reflections, batch_offsets):
  """Assign a 'batch' column to the reflection table"""
  for batch_offset, refl in zip(batch_offsets, reflections):
    xdet, ydet, zdet = [flex.double(x) for x in refl['xyzobs.px.value'].parts()]
    # compute BATCH values - floor() to get (fortran) image captured within
    #                        +1     because FORTRAN counting; zdet+1=image_index
    #                        +off   because            image_index+o=batch
    refl['batch'] = (flex.floor(zdet).iround() + 1) + batch_offset
  return reflections

def get_batch_ranges(experiments, batch_offsets):
  """Get batch ranges for a list of experiments and offsets"""
  batch_ranges = []
  assert len(experiments) == len(batch_offsets)
  image_ranges = get_image_ranges(experiments)
  for batch_offset, image_range in zip(batch_offsets, image_ranges):
    batch_ranges.append((image_range[0] + batch_offset,
      image_range[1] + batch_offset))
  return batch_ranges

def get_image_ranges(experiments):
  """Get image ranges for a list of experiments (including scanless exp.)"""
  # Note, if set to 1,1,for scanless experiments then first batch offset in
  # _calculate_batch_offsets is zero below, bad!
  return [e.scan.get_image_range() if e.scan else (0, 0) for e in experiments]

def calculate_batch_offsets(experiment_list):
  """Take a list of experiments and resolve and return the batch offsets.
  First adds an image_range property as not all experiments have scans."""
  image_ranges = get_image_ranges(experiment_list)
  offsets = _calculate_batch_offsets(image_ranges)
  return offsets

def set_batch_offsets(experiment_list, batch_offsets):
  """Set batch offsets in scan objects. Don't need to set anything for
  scanless experiments, as these are not used with the batch system."""
  for exp, offset in zip(experiment_list, batch_offsets):
    if exp.scan:
      exp.scan.set_batch_offset(offset)

def _calculate_batch_offsets(image_ranges):
  """Take a list of (modified) experiments and resolve and return the batch
  offsets.

  This is the number added to the image number to give the
  batch number, such that:
  - Each experiment has a unique, nonoverlapping, nonconsecutive range
  - None are zero
  - Image number ranges are kept if at all possible
  """

  experiments_to_shift = []
  existing_ranges = set()
  maximum_batch_number = 0
  batch_offsets = [0]*len(image_ranges)

  # Handle zeroth shifts and kept ranges
  for i, image_range in enumerate(image_ranges):
    ilow, ihigh = image_range
    # Check assumptions
    assert ilow <= ihigh, "Inverted image order!?"
    assert ilow >= 0, "Negative image indices are not expected"
    # Don't emit zero: Causes problems with C/fortran number conversion
    if ilow == 0:
      ilow, ihigh = ilow+1, ihigh+1
    # If we overlap with anything, then process later
    if any( ilow <= high+1 and ihigh >= low-1 for low, high in existing_ranges):
      experiments_to_shift.append((i, image_range))
    else:
      batch_offsets[i] = ilow-image_range[0]
      existing_ranges.add((ilow, ihigh))
      maximum_batch_number = max(maximum_batch_number, ihigh)
  # Now handle all the experiments that overlapped by pushing them higher
  for i, image_range in experiments_to_shift:
    start_number = _next_epoch(maximum_batch_number)
    range_width = image_range[1]-image_range[0]+1
    end_number = start_number + range_width - 1
    batch_offsets[i] = start_number - image_range[0]
    maximum_batch_number = end_number
  return batch_offsets

def _next_epoch(val):
  """Find the next number above the existing value that ends in 1, that is
  not consecutive with the current value."""
  if val % 100 == 99:
    return val + 2
  elif val % 100 == 0:
    return val + 101
  else:
    rem = val % 100
    return val - rem + 101
