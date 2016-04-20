from __future__ import division

def run():
  from libtbx import easy_run
  import libtbx.load_env
  import os.path
  from glob import glob

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)

  images = sorted(glob(os.path.join(data_dir, "centroid*.cbf")))
  images_1 = images[0:int(len(images)/2)]
  images_2 = images[int(len(images)/2):]

  easy_run.fully_buffered(
    command=" ".join([
      "dials.import",
      " ".join(images_1),
      "output.datablock=datablock-1.json"])).raise_if_errors()

  result = easy_run.fully_buffered(
    command=" ".join([
      'dials.find_spots',
      'datablock-1.json',
      'output.reflections=strong-1.pickle'])).raise_if_errors()

  easy_run.fully_buffered(
    command=" ".join([
      "dials.import",
      " ".join(images_2),
      "output.datablock=datablock-2.json"])).raise_if_errors()

  result = easy_run.fully_buffered(
    command=" ".join([
      'dials.find_spots',
      'datablock-2.json',
      'output.reflections=strong-2.pickle'])).raise_if_errors()

  result = easy_run.fully_buffered(
    command=" ".join([
      'dials.combine_found_spots',
      'datablock-1.json',
      'datablock-2.json',
      'strong-1.pickle',
      'strong-2.pickle',
      'output.reflections=combined.pickle',
      'output.datablock=combined.json',
      '|'
      'tee',
      'output.txt'])).raise_if_errors()


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print 'OK'
