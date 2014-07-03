from __future__ import division

import os
from iotbx.reflection_file_reader import any_reflection_file
import iotbx.merging_statistics


def run(args):
  for file_name in args:
    # read the file
    reader = any_reflection_file(file_name)

    # now ask for the reflection data as a list of miller arrays
    # we want to get the data unmerged if possible
    as_miller_arrays = reader.as_miller_arrays(merge_equivalents=False)

    # let's see what data we have
    for ma in as_miller_arrays:
      ma.show_summary()

    # look for data arrays with labels I,SIGI
    intensities = None
    for ma in as_miller_arrays:
      if ma.info().labels == ['I', 'SIGI']:
        intensities = ma
        break

    # make sure that the miller array knows it is anomalous and an intensity array
    intensities = intensities.customized_copy(anomalous_flag=True).set_info(
      intensities.info())
    intensities.set_observation_type_xray_intensity()

    # merge the intensities
    merged_intensities = intensities.merge_equivalents().array()
    mtz_dataset = merged_intensities.as_mtz_dataset(column_root_label="I")
    mtz_object = mtz_dataset.mtz_object()
    mtz_object.write(
      os.path.splitext(os.path.basename(file_name))[0] + "_merged.mtz")

    # calculate some merging statistics - this is essentially phenix.merging_statistics
    result = iotbx.merging_statistics.dataset_statistics(intensities)
    result.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
