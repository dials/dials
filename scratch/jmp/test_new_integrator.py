



if __name__ == '__main__':

  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  import sys
  from os.path import join

  # path = '/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data'
  path = '/home/upc86896/Data/Data/i04-BAG-training/dials_processed/'

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
  rlist.compute_zeta_multi(exlist)
  rlist.compute_d(exlist)
  print ""

  from dials.algorithms.integration.interface import Integrator3D
  from dials.algorithms.integration.interface import phil_scope
  from libtbx import phil

  user_phil = phil.parse('''
    mp.max_procs = %d
    block.size=5
    mp.max_tasks = 4
    filter.ice_rings.filter=False
  ''' % nproc)

  params = phil_scope.fetch(source=user_phil).extract()

  integrator = Integrator3D(exlist, rlist, params)
  result = integrator.integrate()
