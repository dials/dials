# LIBTBX_SET_DISPATCHER_NAME dials.background

from __future__ import division

import iotbx.phil

help_message = '''

Examples::

  dials.background image_*.cbf

'''

phil_scope = iotbx.phil.parse("""\
bins = 300
  .type = int
""", process_includes=True)

def main():
  import sys
  run(sys.argv[1:])

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] image_*.cbf" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0 and len(experiments) == 0 and len(reflections) == 0:
    parser.print_help()
    exit()

  assert(len(datablocks) == 1)

  datablock = datablocks[0]
  imagesets = datablock.extract_imagesets()

  assert(len(imagesets) == 1)

  imageset = imagesets[0]

  for indx in imageset.indices():
    print 'For frame %d:' % indx
    background(imageset, indx, params)

def background(imageset, indx, params):
  from dials.array_family import flex
  from libtbx.phil import parse
  from scitbx import matrix
  import math

  detector = imageset.get_detector()
  beam = imageset.get_beam()
  assert(len(detector) == 1)
  detector = detector[0]
  trusted = detector.get_trusted_range()

  n = matrix.col(detector.get_normal()).normalize()
  b = matrix.col(beam.get_s0()).normalize()
  wavelength = beam.get_wavelength()

  if math.fabs(b.dot(n)) < 0.95:
    from libtbx.utils import Sorry
    raise Sorry('Detector not perpendicular to beam')

  data = imageset.get_raw_data(indx)
  assert(len(data) == 1)
  data = data[0]
  negative = (data < 0)
  hot = (data > int(round(trusted[1])))
  bad = negative | hot

  from dials.algorithms.spot_finding.factory import SpotFinderFactory
  from dials.algorithms.spot_finding.factory import phil_scope

  data = data.as_double()

  spot_params = phil_scope.fetch(source=parse("")).extract()
  threshold_function = SpotFinderFactory.configure_threshold(spot_params)
  peak_pixels = threshold_function.compute_threshold(data, ~bad)
  signal = data.select(peak_pixels.iselection())
  background = data.select((~bad & ~peak_pixels).iselection())

  # print some summary information
  print 'Mean background: %.3f' % (flex.sum(background) / background.size())
  print 'Max/total signal pixels: %.0f / %.0f' % (flex.max(signal),
                                                 flex.sum(signal))
  print 'Peak/background/hot pixels: %d / %d / %d' % (peak_pixels.count(True),
                                                      background.size(),
                                                      hot.count(True))

  # compute histogram of two-theta values, then same weighted
  # by pixel values, finally divide latter by former to get
  # the radial profile out, need to set the number of bins
  # sensibly; flex.histogram does not allow weights so use
  # numpy.histogram to get the same effect... inspired by
  # method in PyFAI

  two_theta_array = detector.get_two_theta_array(beam.get_s0())
  two_theta_array.reshape(flex.grid(data.focus()))
  two_theta_array.as_1d().set_selected((bad | peak_pixels).iselection(), 0.0)
  data.as_1d().set_selected((bad | peak_pixels).iselection(), 0.0)

  # numpy land :-(
  import numpy
  n_bins = params.bins
  n_two_theta = two_theta_array.as_1d().as_numpy_array()
  n_data = data.as_1d().as_numpy_array()
  n_data2 = (data * data).as_1d().as_numpy_array()

  tt_range = 0, flex.max(two_theta_array)

  ref, junk = numpy.histogram(n_two_theta, bins=n_bins, range=tt_range)
  val, junk = numpy.histogram(n_two_theta, bins=n_bins, weights=n_data,
                              range=tt_range)
  val2, junk = numpy.histogram(n_two_theta, bins=n_bins, weights=n_data2,
                               range=tt_range)

  I = val / ref
  I2 = val2 / ref

  print '%8s %8s %8s' % ('d', 'I', 'sig')
  for j in range(0, len(I)):
    tt = tt_range[0] + (j + 0.5) * (1.0 / n_bins) * (tt_range[1] - tt_range[0])
    d = wavelength / (2.0 * math.sin(0.5 * tt))
    print '%8.3f %8.3f %8.3f' % (d, I[j], math.sqrt(I2[j] - I[j]**2))

if __name__ == '__main__':
  main()
