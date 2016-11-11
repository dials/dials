# LIBTBX_SET_DISPATCHER_NAME dev.dials.image_correlation

from __future__ import division

import iotbx.phil
from dials.array_family import flex

help_message = '''

Examples::

  dev.dials.image_correlation datablock.json image=1,2,3

'''

phil_scope = iotbx.phil.parse("""\
  image=None
    .type = ints
    .help = "image numbers to analyse"
""", process_includes=True)

def main():
  import sys
  run(sys.argv[1:])

def extract_signal_mask(data):
  from dials.algorithms.spot_finding.factory import SpotFinderFactory
  from dials.algorithms.spot_finding.factory import phil_scope

  data = data.as_double()

  from dxtbx import datablock

  spot_params = phil_scope.fetch(source=iotbx.phil.parse(
    "spotfinder.threshold.xds.gain=1")).extract()
  threshold_function = SpotFinderFactory.configure_threshold(
    spot_params, None)
  negative = (data < 0)
  signal = threshold_function.compute_threshold(data, ~negative)

  return signal

def image_correlation(a, b):

  sig_a = extract_signal_mask(a)
  sig_b = extract_signal_mask(b)

  a = a.as_1d().as_double()
  b = b.as_1d().as_double()
  sel = (a > 0) & (b > 0) & sig_a & sig_b
  _a = a.select(sel)
  _b = b.select(sel)
  return sel.count(True), flex.linear_correlation(x=a, y=b).coefficient()

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks
  import libtbx.load_env

  usage = "%s [options] datablock.json" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    exit()

  assert(len(datablocks) == 1)

  datablock = datablocks[0]
  imagesets = datablock.extract_imagesets()

  assert(len(imagesets) == 1)

  imageset = imagesets[0]

  images = params.image

  for j, img_a in enumerate(images[:-1]):
    for img_b in images[j+1:]:
      a = imageset.get_raw_data(img_a)[0]
      b = imageset.get_raw_data(img_b)[0]
      n, cc = image_correlation(a, b)
      print '%5d %5d %7d %.4f' % (img_a, img_b, n, cc)

if __name__ == '__main__':
  main()
