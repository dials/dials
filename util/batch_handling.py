"""
Functions to help with calculating batch properties for experiments objects.
"""
from math import ceil, floor, log
import logging
import copy as copy
from dials.array_family import flex

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

def are_all_batch_ranges_set_by_scaling(experiments):
  for experiment in experiments:
    if experiment.scaling_model:
      if not 'valid_batch_range' in experiment.scaling_model.configdict:
        return False
    else:
      return False
  return True

def calculate_batch_offsets(experiment_list):
  """Take a list of experiments and resolve and return the batch offsets.
  First adds an image_range property as not all experiments have scans."""
  if are_all_batch_ranges_set_by_scaling(experiment_list):
    return [e.scan.get_batch_offset() for e in experiment_list]
  image_ranges = get_image_ranges(experiment_list)
  offsets = _calculate_batch_offsets(image_ranges)
  return offsets

def set_batch_offsets(experiment_list, batch_offsets):
  """Set batch offsets in scan objects. Don't need to set anything for
  scanless experiments, as these are not used with the batch system."""
  for exp, offset in zip(experiment_list, batch_offsets):
    if exp.scan:
      exp.scan.set_batch_offset(offset)

def exclude_batches_in_reflections(reflections, experiments,
  in_use_batch_ranges, exclude_batches):
  """Determine the new batch ranges when some batches are excluded.
  Exclude batches is a list of tuples e.g. [(101,200), (301,400)].
  Also set the relevant flag in the reflection tables."""
  valid_batch_ranges = [] #list of valid batch ranges after excluding
  n_ranges_excluded = 0
  for exp, refl, batch_range in zip(experiments, reflections, in_use_batch_ranges):
    excluded_from_dataset_i = []
    for exclude_pair in exclude_batches:
      if list(batch_range) == exclude_pair: #raise, else very tricky to implement!
        raise ValueError("Trying to exclude a whole dataset")
      if exclude_pair[1] <= batch_range[1] and exclude_pair[0] >= batch_range[0]:
        mask1 = refl['batch'] <= exclude_pair[1]
        mask2 = refl['batch'] >= exclude_pair[0]
        refl.set_flags(mask1 & mask2, refl.flags.user_excluded_in_scaling)
        logger.info("""Excluded batch range %s for dataset %s""",
          exclude_pair, exp.identifier)
        excluded_from_dataset_i.extend(list(exclude_pair))
        n_ranges_excluded += 1
    if excluded_from_dataset_i:
      #calculate min-max range that valid batches span (may include 'gaps')
      new_batch_0, new_batch_1 = batch_range
      s = sorted(excluded_from_dataset_i)
      if batch_range[0] == s[0]:
        new_batch_0 = s[1]+1
      if batch_range[1] == s[-1]:
        new_batch_1 = s[-2]-1
      valid_batch_ranges.append((new_batch_0, new_batch_1))
    else:
      valid_batch_ranges.append(batch_range)
  if n_ranges_excluded != len(exclude_batches):
    logger.warn("Warning: not all requested ranges were excluded, please check input")
  return reflections, valid_batch_ranges

def get_current_batch_ranges_for_scaling(experiments, batch_ranges):
  """Inspect scaling model to see if the batch range has already been reduced"""
  in_use_batch_ranges = copy.deepcopy(batch_ranges)
  for i, experiment in enumerate(experiments):
    if experiment.scaling_model is not None:
      if 'valid_batch_range' in experiment.scaling_model.configdict:
        in_use_batch_ranges[i] = experiment.scaling_model.configdict['valid_batch_range']
  return in_use_batch_ranges

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
