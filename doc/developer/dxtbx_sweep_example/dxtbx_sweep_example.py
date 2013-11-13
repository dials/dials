from __future__ import division


if __name__ == '__main__':

  import os
  import libtbx.load_env
  from glob import glob
  from dxtbx.imageset import ImageSetFactory

  # Check dials_regression is configured
  try:
    path = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    raise

  # Find the filenames
  template = os.path.join(path, 'centroid_test_data', 'centroid_*.cbf')
  filenames = glob(template)

  # Create the sweep
  sweep = ImageSetFactory.new(filenames)
  assert(len(sweep) == 1)
  sweep = sweep[0]

  # Get the models
  beam = sweep.get_beam()
  detector = sweep.get_detector()
  gonio = sweep.get_goniometer()
  scan = sweep.get_scan()
  print beam
  print detector
  print gonio
  print scan

  print "sweep: ", sweep
  print "sweep indices: ", sweep.indices()
  print "sweep array range: ", sweep.get_array_range()

  # Get a sub sweep
  sub_sweep = sweep[3:6]
  print "sub_sweep: ", sub_sweep
  print "sub_sweep indices: ", sub_sweep.indices()
  print "sub_sweep array range: ", sub_sweep.get_array_range()

  # Loop through sub sweep
  offset = sub_sweep.get_array_range()[0]
  for i, f in enumerate(sub_sweep):
    print "Image {0} has size {1}".format(i + offset, f.all())

  # Extract a volume
  volume = sweep.to_array()
  print "sweep volume 1 size: ", volume.all()

  volume = sweep.to_array((2, 7))
  print "sweep volume 2 size: ", volume.all()

  volume = sweep[2:7].to_array()
  print "sweep volume 3 size: ", volume.all()

  volume = sweep.to_array((2, 7, 100, 200, 100, 300))
  print "sweep volume 4 size: ", volume.all()

  print "Format: ", sweep.reader().get_format()
  print "Template: ", sweep.get_template()
