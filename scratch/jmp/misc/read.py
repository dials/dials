


if __name__ == '__main__':

  from dials.array_family import flex
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory
  import os.path
  path = "/home/upc86896/Projects/cctbx/sources/dials_regression/centroid_test_data"

  rlist_filename = os.path.join(path, "integrated.pickle")
  exlist_filename = os.path.join(path, "experiments.json")

  rlist = flex.reflection_table.from_pickle(rlist_filename)
  exlist = ExperimentListFactory.from_json_file(exlist_filename)

  panel = rlist['panel']
  bbox = rlist['bbox']

  rlist['shoebox'] = flex.shoebox(panel, bbox)
  rlist['shoebox'].allocate()

  rlist.fill_shoeboxes(exlist[0].imageset)
