"""
Functions to help with calculating batch properties for experiments objects.
"""
from math import ceil, floor, log
from dials.array_family import flex

def assign_batches_to_reflections(reflections, batch_offsets):

  for batch_offset, refl in zip(batch_offsets, reflections):
    xdet, ydet, zdet = [flex.double(x) for x in refl['xyzobs.px.value'].parts()]
    # compute BATCH values - floor() to get (fortran) image captured within
    #                        +1     because FORTRAN counting; zdet+1=image_index
    #                        +off   because            image_index+o=batch
    refl['batch'] = (flex.floor(zdet).iround() + 1) + batch_offset
  return reflections

def get_batch_ranges(experiments, batch_offsets):
  batch_ranges = []
  assert len(experiments) == len(batch_offsets)
  for i, experiment in enumerate(experiments):
    assign_image_range_to_experiment(experiment)
    batch_ranges.append((experiment.image_range[0] + batch_offsets[i],
      experiment.image_range[1] + batch_offsets[i]))
  return batch_ranges

def get_image_ranges(experiments):
  #experiments must be passed in as list(experiments), not an ExperimentList
  image_ranges = []
  for experiment in experiments:
    assign_image_range_to_experiment(experiment)
    image_ranges.append(experiment.image_range)
  return image_ranges

def calculate_batch_offsets(experiment_list):
  """Take a list of experiments and resolve and return the batch offsets.
  First adds an image_range property as not all experiments have scans."""

  for experiment in experiment_list:
    assign_image_range_to_experiment(experiment)
  offsets = _calculate_batch_offsets(experiment_list)
  return offsets

def set_batch_offsets(experiment_list, batch_offsets):
  """Set batch offsets in scan objects. Don't need to set anything for
  scanless experiments, as these are not used with the batch system."""
  for exp, offset in zip(experiment_list, batch_offsets):
    if exp.scan:
      exp.scan.set_batch_offset(offset)

def assign_image_range_to_experiment(experiment):
  if not hasattr(experiment, 'image_range'):
    if experiment.scan:
      experiment.image_range = experiment.scan.get_image_range()
    else:
      experiment.image_range = 0, 0
      # Note, if set to 1,1, then first batch offset is zero below, bad!

def calculate_new_batch_ranges(in_use_batch_ranges, exclude_batches):
  """Determine the new batch ranges.
  Exclude batches is a list of tuples e.g. [(101,200), (301,400)]"""
  valid_batch_ranges = [] #list of valid batch ranges after excluding
  n_ranges_excluded = 0
  for batch_range in in_use_batch_ranges:
    excluded_from_dataset_i = []
    for exclude_pair in exclude_batches:
      if list(batch_range) == exclude_pair: #raise, else very tricky to implement!
        raise ValueError("Trying to exclude a whole dataset")
      if exclude_pair[1] <= batch_range[1] and exclude_pair[0] >= batch_range[0]:
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
  return valid_batch_ranges

def exclude_batches_in_reflections(reflections, experiments,
  valid_batch_ranges, in_use_batch_ranges):
  """Exclude batches from the tables by setting user_excluded_in_scaling flag."""
  for (exp, refl, valid, in_use) in zip(experiments, reflections,
    valid_batch_ranges, in_use_batch_ranges):
    if valid != in_use:
      assert valid[0] >= in_use[0]
      assert valid[1] <= in_use[1]
      mask1 = refl['batch'] < valid[0]
      mask2 = refl['batch'] > valid[1]
      refl.set_flags(mask1 & mask2, refl.flags.user_excluded_in_scaling)
      logger.info("""Reduced batch range for dataset %s to %s""",
        exp.identifier, valid)
  return reflections

def get_current_batch_ranges_for_scaling(experiments, batch_ranges):
  """Inspect scaling model to see if the batch range has already been reduced"""
  in_use_batch_ranges = copy.deepcopy(batch_ranges)
  for i, experiment in enumerate(experiments):
    if experiment.scaling_model is not None:
      if 'valid_batch_range' in experiment.scaling_model.configdict:
        in_use_batch_ranges[i] = experiment.scaling_model.configdict['valid_batch_range']
  return in_use_batch_ranges

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
      experiment.batch_offset = ilow-experiment.image_range[0]
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
  """Find the next number above the existing value that ends in 1, that is
  not consecutive with the current value."""
  if val % 100 == 99:
    return val + 2
  elif val % 100 == 0:
    return val + 101
  else:
    rem = val % 100
    return val - rem + 101
