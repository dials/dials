
from __future__ import division


if __name__ == '__main__':

  from dials.util.options import OptionParser
  from dials.algorithms.profile_model.profile_model import ProfileModelList, phil_scope
  from libtbx.phil import parse
  import sys
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  from dials.array_family import flex
  from dials.viewer.tools import show_reflection

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
