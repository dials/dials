


if __name__ == '__main__':

  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  import sys

  assert(len(sys.argv) > 1)
  filename = sys.argv[1]
  if len(sys.argv) > 2:
    num_proc = int(sys.argv[2])
  else:
    num_proc = 1

  from math import pi
  exlist = ExperimentListFactory.from_json_file(sys.argv[1])
  rlist = flex.reflection_table.from_predictions(exlist[0])
  rlist['id'] = flex.int(len(rlist), 0)
  rlist.compute_bbox(exlist[0], nsigma=3, sigma_d=0.024*pi/180,
                     sigma_m=0.044*pi/180)

  from dials.algorithms.integration.fast_integrator import integrate_quickly

  integrate_quickly(exlist, rlist, num_proc=num_proc)
