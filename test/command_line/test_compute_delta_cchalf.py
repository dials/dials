
from __future__ import absolute_import, division, print_function
from iotbx.reflection_file_reader import any_reflection_file
from dials.algorithms.statistics.delta_cchalf import PerImageCChalfStatistics
from os.path import join

def test_compute_delta_cchalf(dials_regression):

  filename = join(dials_regression, "delta_cchalf_test_data", "test.XDS_ASCII.mtz")

  # Read the mtz file
  reader = any_reflection_file(filename)

  # Get the columns as miller arrays
  miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

  # Select the desired columns
  intensities = None
  batches = None
  for array in miller_arrays:
    if array.info().labels == ['I', 'SIGI']:
      intensities = array
    if array.info().labels == ['BATCH']:
      batches = array
  assert intensities is not None
  assert batches is not None
  assert len(batches.data()) == len(intensities.data())

  # Get the unit cell and space group
  unit_cell = intensities.unit_cell()
  space_group = intensities.crystal_symmetry().space_group()

  # The reflection data
  miller_index = intensities.indices()
  batch = batches.data()
  intensity = intensities.data()
  variance = intensities.sigmas()**2

  # Create unit cell list
  min_batch = min(batch)
  dataset = batch - min_batch
  num_datasets = max(dataset)+1
  unit_cell_list = [unit_cell for i in range(num_datasets)]

  # Compute the CC 1/2 Stats
  statistics = PerImageCChalfStatistics(
    miller_index,
    dataset,
    intensity,
    variance,
    unit_cell_list,
    space_group,
    nbins=1)

  mean_cchalf = statistics.mean_cchalf()
  cchalf_i = statistics.cchalf_i()

  assert abs(100*mean_cchalf - 94.582) < 1e-3
  assert abs(100*cchalf_i[0] - 79.587) < 1e-3
  assert abs(100*cchalf_i[1] - 94.238) < 1e-3
