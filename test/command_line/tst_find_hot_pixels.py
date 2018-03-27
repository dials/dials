from __future__ import absolute_import, division
import os
from dials.array_family import flex # import dependency
import cPickle as pickle
import libtbx.load_env
from libtbx import easy_run
from glob import glob

class Test(object):

  def __init__(self):

    # Check we have dials_regression
    if not libtbx.env.has_module("dials_regression"):
      print "Skipping: dials_regression not present"
      exit(0)

    # set the directory
    self.directory = libtbx.env.find_in_repositories(
      relative_path="dials_regression/centroid_test_data",
      test=os.path.isdir)

  def run(self):

    from os.path import join, exists

    template = glob(join(self.directory, "centroid*.cbf"))
    args = [
      "dials.find_spots", ' '.join(template),
      "output.datablock=datablock.json",
      "output.reflections=spotfinder.pickle",
      "output.shoeboxes=True"
    ]
    result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()

    assert exists("datablock.json")
    assert exists("spotfinder.pickle")

    args = [
      "dials.find_hot_pixels",
      "input.datablock=datablock.json",
      "input.reflections=spotfinder.pickle",
      "output.mask=hot_mask.pickle"
    ]
    result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()

    assert exists("hot_mask.pickle")


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
