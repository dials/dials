



if __name__ == '__main__':

  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.algorithms.profile_model.profile_model import ProfileModelList
  from dials.algorithms.profile_model.profile_model import ProfileModel
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

  print "Reading Experiments"
  from math import pi
  experiments = ExperimentListFactory.from_json_file(experiment_list_filename)

  profile_model = ProfileModelList()
  profile_model.append(ProfileModel(
    n_sigma=3,
    sigma_b=0.024*pi/180,
    sigma_m=0.044*pi/180))


  print "Predicting Reflections"
  rlist = flex.reflection_table.from_predictions(experiments[0])
  rlist['id'] = flex.int(len(rlist), 0)
  rlist.compute_bbox(experiments, profile_model)
  rlist.compute_zeta_multi(experiments)
  rlist.compute_d(experiments)
  print ""

  print "Creating params"
  from dials.algorithms.integration.integrator import IntegratorFactory
  from dials.algorithms.integration.integrator import phil_scope
  from libtbx import phil

  user_phil = phil.parse('''
    integration {
      mp.max_procs = %d
      block.size=5
      filter.ice_rings.filter=False
      intensity.algorithm=sum3d
    }
  ''' % nproc)

  working_phil = phil_scope.fetch(source=user_phil)
  params = working_phil.extract()

  print "Integrating"
  integrator = IntegratorFactory.create(params, experiments, profile_model, rlist)
  result = integrator.integrate()

  result.as_pickle("temp.pickle")
