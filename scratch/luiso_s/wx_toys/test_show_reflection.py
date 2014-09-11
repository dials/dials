
from __future__ import division


if __name__ == '__main__':

  from dials.util.options import OptionParser
  from dials.algorithms.profile_model.profile_model import ProfileModelList, phil_scope
  from libtbx.phil import parse
  import sys
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  from dials.viewer.tools import show_reflection

  from os.path import join
  import libtbx.load_env
  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    exit(0)
  path = join(dials_regression, "centroid_test_data")
  import sys
  assert(len(sys.argv) == 1)
  sys.argv.append(join(path, "experiments.json"))
  sys.argv.append(join(path, "profile.phil"))

  parser = OptionParser(phil=phil_scope)
  params, options, args = parser.parse_args()
  assert(len(args) == 1)

  exlist = ExperimentListFactory.from_json_file(args[0])
  assert(len(exlist) == 1)


  profile_model = ProfileModelList.load(params)

  rlist = flex.reflection_table.from_predictions_multi(exlist)
  rlist.compute_bbox(exlist, profile_model)
  rlist['shoebox'] = flex.shoebox(rlist['panel'], rlist['bbox'])
  rlist['shoebox'].allocate()

  rlist.extract_shoeboxes(exlist[0].imageset)

  show_reflection(rlist[len(rlist)//2])
  #show_reflection(rlist[len(rlist)//2], orient = "porTrait")
  #show_reflection(rlist[len(rlist)//2], orient = "lanDscape")
