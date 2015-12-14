#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2015) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

"""
Development testing of restraints_parameterisation and associated classes

"""

# Python and cctbx imports
from __future__ import division
import os
from libtbx.phil import parse
import libtbx.load_env # required for libtbx.env.find_in_repositories
from dxtbx.model.experiment.experiment_list import ExperimentListFactory
from dials.algorithms.refinement import RefinerFactory
from dials.array_family import flex

# The phil scope
from dials.algorithms.refinement.refiner import phil_scope
user_phil = parse('''
refinement
{
  parameterisation
  {
    crystal
    {
      unit_cell
      {
        restraints
        {
          tie_to_target
          {
            values=10,10,10,90,90,90
            sigmas=1,1,1,1,1,1
            id=0
          }
        }
      }
    }
  }
}
''')

working_phil = phil_scope.fetch(source=user_phil)
working_params = working_phil.extract()

def test1(params):

  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)

  # use the multi stills test data
  data_dir = os.path.join(dials_regression, "refinement_test_data", "multi_stills")
  experiments_path = os.path.join(data_dir, "combined_experiments.json")
  pickle_path = os.path.join(data_dir, "combined_reflections.pickle")

  experiments = ExperimentListFactory.from_json_file(experiments_path,
                check_format=False)
  reflections = flex.reflection_table.from_pickle(pickle_path)

  refiner = RefinerFactory.from_parameters_data_experiments(params,
        reflections, experiments)

  # hack to extract the restraints parameterisation from the Refiner
  rp = refiner._target._restraints_parameterisation

  # enter interactive console
  from dials.util.command_line import interactive_console; interactive_console()

if __name__ == '__main__':

  test1(working_params)


