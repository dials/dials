



if __name__ == '__main__':

  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  import sys
  from os.path import join

  path = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'

  experiment_list_filename = join(path, "experiments.json")

  if len(sys.argv) > 1:
    nproc = int(sys.argv[1])
  else:
    nproc = 1

  from math import pi
  exlist = ExperimentListFactory.from_json_file(experiment_list_filename)
  rlist = flex.reflection_table.from_predictions(exlist[0])
  rlist['id'] = flex.size_t(len(rlist), 0)
  rlist.compute_bbox(exlist[0], nsigma=3, sigma_d=0.024*pi/180,
                     sigma_m=0.044*pi/180)

  from dials.algorithms.integration.interface import Integrator3D

  integrator = Integrator3D(exlist, rlist, num_tasks=nproc, max_procs=5)
  result = integrator.integrate()
