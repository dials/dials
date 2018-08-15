# LIBTBX_SET_DISPATCHER_NAME dials.background
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division, print_function

import iotbx.phil
from scitbx.array_family import flex
from libtbx.utils import Sorry

help_message = '''

Examples::

  dials.background image_*.cbf

'''

phil_scope = iotbx.phil.parse("""\
n_bins = 100
  .type = int
images = None
  .type = ints
  .help = "Images on which to perform the analysis (otherwise use all images)"
plot = False
  .type = bool

masking {
  include scope dials.util.masking.phil_scope
}
""", process_includes=True)

def main():
  import sys
  run(sys.argv[1:])

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks, flatten_experiments
  import libtbx.load_env

  usage = "%s [options] image_*.cbf" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    read_experiments=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  # Ensure we have either a data block or an experiment list
  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)
  if [len(datablocks), len(experiments)].count(1) != 1:
      self.parser.print_help()
      print("Please pass either a datablock or an experiment list\n")
      return

  if datablocks:
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
  else:
    imagesets = experiments.imagesets()

  if len(imagesets) != 1:
    raise Sorry("Please pass a datablock or experiment list that contains a "
      "single imageset")
  imageset = imagesets[0]

  first, last = imageset.get_scan().get_image_range()
  images = range(first, last + 1)

  if params.images:
    if min(params.images) < first or max(params.images) > last:
      raise Sorry("image outside of scan range")
    images = params.images

  d_spacings = []
  intensities = []
  sigmas = []

  for indx in images:
    print('For image %d:' % indx)
    indx -= first # indices passed to imageset.get_raw_data start from zero
    d, I, sig = background(imageset, indx, n_bins=params.n_bins,
        mask_params=params.masking)

    print('%8s %8s %8s' % ('d', 'I', 'sig'))
    for j in range(len(I)):
      print('%8.3f %8.3f %8.3f' % (d[j], I[j], sig[j]))

    d_spacings.append(d)
    intensities.append(I)
    sigmas.append(sig)

  if params.plot:
    from matplotlib import pyplot
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.set_xlabel(r'resolution ($\AA$)')
    ax.set_ylabel(r'$\langle I_b \rangle$')
    for d, I, sig in zip(d_spacings, intensities, sigmas):
      ds2 = 1/flex.pow2(d)
      ax.plot(ds2, I)
    xticks = ax.get_xticks()
    from math import sqrt
    x_tick_labs = ["" if e <= 0.0 else "{:.2f}".format(sqrt(1./e))
                   for e in xticks]
    ax.set_xticklabels(x_tick_labs)

    pyplot.show()

def background(imageset, indx, n_bins, mask_params=None):
  from dials.array_family import flex
  from libtbx.phil import parse
  from scitbx import matrix
  import math

  if mask_params is None:
    # Default mask params for trusted range
    mask_params = phil_scope.fetch(parse("")).extract().masking

  from dials.util.masking import MaskGenerator
  mask_generator = MaskGenerator(mask_params)
  mask = mask_generator.generate(imageset)

  detector = imageset.get_detector()
  beam = imageset.get_beam()
  # Only working with single panel detector for now
  assert(len(detector) == 1)
  panel = detector[0]
  mask = mask[0]

  n = matrix.col(panel.get_normal()).normalize()
  b = matrix.col(beam.get_s0()).normalize()
  wavelength = beam.get_wavelength()

  if math.fabs(b.dot(n)) < 0.95:
    raise Sorry('Detector not perpendicular to beam')

  data = imageset.get_raw_data(indx)
  assert(len(data) == 1)
  data = data[0]

  from dials.algorithms.spot_finding.factory import SpotFinderFactory
  from dials.algorithms.spot_finding.factory import phil_scope as spot_phil

  data = data.as_double()

  from dxtbx import datablock

  spot_params = spot_phil.fetch(source=parse("")).extract()
  threshold_function = SpotFinderFactory.configure_threshold(
    spot_params, datablock.DataBlock([imageset]))
  peak_pixels = threshold_function.compute_threshold(data, mask)
  signal = data.select(peak_pixels.iselection())
  background_pixels = (mask & ~peak_pixels)
  background = data.select(background_pixels.iselection())

  # print some summary information
  print('Mean background: %.3f' % (flex.sum(background) / background.size()))
  if len(signal) > 0:
    print('Max/total signal pixels: %.0f / %.0f' % (flex.max(signal),
                                                    flex.sum(signal)))
  else:
    print('No signal pixels on this image')
  print('Peak/background/masked pixels: %d / %d / %d' % (peak_pixels.count(True),
                                                      background.size(),
                                                      mask.count(False)))

  # compute histogram of two-theta values, then same weighted
  # by pixel values, finally divide latter by former to get
  # the radial profile out, need to set the number of bins
  # sensibly; inspired by method in PyFAI

  two_theta_array = panel.get_two_theta_array(beam.get_s0())
  two_theta_array = two_theta_array.as_1d().select(background_pixels.iselection())

  # Use flex.weighted_histogram
  h0 = flex.weighted_histogram(two_theta_array, n_slots=n_bins)
  h1 = flex.weighted_histogram(two_theta_array, background, n_slots=n_bins)
  h2 = flex.weighted_histogram(two_theta_array, background * background, n_slots=n_bins)

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
