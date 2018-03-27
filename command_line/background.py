# LIBTBX_SET_DISPATCHER_NAME dials.background
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function

import iotbx.phil
from scitbx.array_family import flex

help_message = '''

Examples::

  dials.background image_*.cbf

'''

phil_scope = iotbx.phil.parse("""\
n_bins = 100
  .type = int
frames = None
  .type = int
  .multiple = True
plot = False
  .type = bool
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

  images = imageset.indices()
  if params.frames:
    images = params.frames

  d_spacings = []
  intensities = []
  sigmas = []

  for indx in images:
    print('For frame %d:' % indx)
    d, I, sig = background(imageset, indx, n_bins=params.n_bins)

    print('%8s %8s %8s' % ('d', 'I', 'sig'))
    for j in range(len(I)):
      print('%8.3f %8.3f %8.3f' % (d[j], I[j], sig[j]))

    d_spacings.append(d)
    intensities.append(I)
    sigmas.append(sig)

  if params.plot:
    from matplotlib import pyplot
    fig = pyplot.figure()
    for d, I, sig in zip(d_spacings, intensities, sigmas):
      ds2 = 1/flex.pow2(d)
      pyplot.plot(ds2, I)

    pyplot.show()

def background(imageset, indx, n_bins):
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

  from dxtbx import datablock

  spot_params = phil_scope.fetch(source=parse("")).extract()
  threshold_function = SpotFinderFactory.configure_threshold(
    spot_params, datablock.DataBlock([imageset]))
  peak_pixels = threshold_function.compute_threshold(data, ~bad)
  signal = data.select(peak_pixels.iselection())
  background = data.select((~bad & ~peak_pixels).iselection())

  # print some summary information
  print('Mean background: %.3f' % (flex.sum(background) / background.size()))
  print('Max/total signal pixels: %.0f / %.0f' % (flex.max(signal),
                                                  flex.sum(signal)))
  print('Peak/background/hot pixels: %d / %d / %d' % (peak_pixels.count(True),
                                                      background.size(),
                                                      hot.count(True)))

  # compute histogram of two-theta values, then same weighted
  # by pixel values, finally divide latter by former to get
  # the radial profile out, need to set the number of bins
  # sensibly; flex.histogram does not allow weights so use
  # numpy.histogram to get the same effect... inspired by
  # method in PyFAI

  data = data.as_1d()
  two_theta_array = detector.get_two_theta_array(beam.get_s0())
  two_theta_array.set_selected((bad | peak_pixels).iselection(), 0.0)
  data.set_selected((bad | peak_pixels).iselection(), 0.0)

  # new fangled flex.weighted_histogram :-)
  h0 = flex.weighted_histogram(two_theta_array, n_slots=n_bins)
  h1 = flex.weighted_histogram(two_theta_array, data, n_slots=n_bins)
  h2 = flex.weighted_histogram(two_theta_array, data * data, n_slots=n_bins)

  d0 = h0.slots()
  d1 = h1.slots()
  d2 = h2.slots()

  I = d1 / d0
  I2 = d2 / d0
  sig = flex.sqrt(I2 - flex.pow2(I))

  tt = h0.slot_centers()
  d_spacings = wavelength / (2.0 * flex.sin(0.5 * tt))

  return d_spacings, I, sig

if __name__ == '__main__':
  main()
