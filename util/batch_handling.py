"""
Functions to help with calculating batch properties for experiments objects.
"""
from math import ceil, floor, log

def calculate_batch_offsets(experiment_list):
  """Take a list of experiments and resolve and return the batch offsets.
  First adds an image_range property as not all experiments have scans."""

  for experiment in experiment_list:
    assign_image_range_to_experiment(experiment)
  offsets = _calculate_batch_offsets(experiment_list)
  return offsets

def assign_image_range_to_experiment(experiment):
  if not hasattr(experiment, 'image_range'):
    if experiment.scan:
      experiment.image_range = experiment.scan.get_image_range()
    else:
      experiment.image_range = 0, 0
      # Note, if set to 1,1, then first batch offset is zero below, bad!

def _calculate_batch_offsets(experiments):
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
  batch_offsets = [0]*len(experiments)

  # Handle zeroth shifts and kept ranges
  for i, experiment in enumerate(experiments):
    ilow, ihigh = experiment.image_range
    # Check assumptions
    assert ilow <= ihigh, "Inverted image order!?"
    assert ilow >= 0, "Negative image indices are not expected"
    # Don't emit zero: Causes problems with C/fortran number conversion
    if ilow == 0:
      ilow, ihigh = ilow+1, ihigh+1
    # If we overlap with anything, then process later
    if any( ilow <= high+1 and ihigh >= low-1 for low, high in existing_ranges):
      experiments_to_shift.append((i, experiment))
    else:
      batch_offsets[i] = ilow-experiment.image_range[0]
      existing_ranges.add((ilow, ihigh))
      maximum_batch_number = max(maximum_batch_number, ihigh)
  # Now handle all the experiments that overlapped by pushing them higher
  for i, experiment in experiments_to_shift:
    start_number = _next_epoch(maximum_batch_number)
    range_width = experiment.image_range[1]-experiment.image_range[0]+1
    end_number = start_number + range_width - 1
    batch_offsets[i] = start_number - experiment.image_range[0]
    maximum_batch_number = end_number
    experiment.batch_offset = batch_offsets[i]
  return batch_offsets

def _next_epoch(val):
  """Find the next number divisible by 10 above existing value, that is
  not consecutive with the current value."""
  if val % 10 == 9:
    return val + 11
  elif val % 10 == 0:
    return val + 10
  else:
    rem = val % 10
    return val - rem + 10
