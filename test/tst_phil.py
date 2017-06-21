from __future__ import absolute_import, division

import os
from dials.util.phil import parse

from mock import Mock, patch

# Modules the phil parser uses, that we want to mock
import dials.array_family.flex
import dxtbx.model.experiment_list
import dxtbx.datablock

class Test(object):
  def __init__(self):
    # Only use these filenames for verification
    self.path = "centroid_test_data"
    self.datablock_path = os.path.join(self.path, "datablock.json")
    self.experiments_path = os.path.join(self.path, "experiments.json")
    self.reflections_path = os.path.join(self.path, "integrated.pickle")

  @patch("os.path.exists", Mock(return_value=True))
  @patch("dials.array_family.flex")
  @patch("dxtbx.model.experiment_list.ExperimentListFactory")
  @patch("dxtbx.datablock.DataBlockFactory")
  def run(self, DataBlockFactory, ExperimentListFactory, flex, *args):

    phil_scope = parse('''
      input {
        reflections = %s
          .type = reflection_table
        datablock = %s
          .type = datablock
        experiments = %s
          .type = experiment_list
      }
    ''' % (self.reflections_path,
           self.datablock_path,
           self.experiments_path))

    params = phil_scope.extract()
    # Check the right filenames were parsed
    assert(params.input.reflections.filename == self.reflections_path)
    assert(params.input.datablock.filename == self.datablock_path)
    assert(params.input.experiments.filename == self.experiments_path)
    # Check that we got the expected objects back
    assert isinstance(params.input.reflections.data, Mock)
    assert isinstance(params.input.datablock.data, Mock)
    assert isinstance(params.input.experiments.data, Mock)
    # Check we had the correct calls made
    flex.reflection_table.from_pickle.assert_called_once_with(self.reflections_path)
    assert DataBlockFactory.from_json_file.call_args[0] == (self.datablock_path,)
    assert ExperimentListFactory.from_json_file.call_args[0] == (self.experiments_path,)

    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
