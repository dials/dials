from __future__ import division

from libtbx.phil import command_line
import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments

import logging
logger = logging.getLogger(__name__)

phil_scope = iotbx.phil.parse("""\
output {
  mask = hot_pixels.pickle
    .type = path
}

verbosity = 1
  .type = int(value_min=0)
  .help = "The verbosity level"
""", process_includes=True)

help_message = '''

  This program looks through the output of dials.find_spots to determine if any
  of the pixels are possible "hot" pixels. If a pixel is "hot" this will mean
  that it will have a consistently high value throughout the dataset. This
  program simply selects all pixels which are labelled as strong on all images
  in the dataset as "hot".

  Note that if you have still data or a small dataset, this is likely to produce
  lots of false positives; however, if you have a large rotation dataset, it is
  likely to be reasonably accurate.

  The program returns a file names hot_pixels.pickle which contains a boolean mask
  with True pixels being OK and False pixels being "hot" pixels.

  Examples::
    dials.find_hot_pixels datablock.json strong.pickle

'''

def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util import log
  import cPickle as pickle
  usage = "%s [options] datablock.json strong.pickle" % \
    libtbx.env.dispatcher_name

  # Create the option parser
  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=True,
    check_format=False,
    epilog=help_message)

  # Get the parameters
  params, options = parser.parse_args(show_diff_phil=False)

  # Configure the log
  log.config(
    params.verbosity,
    info='dials.find_hot_pixels.log',
    debug='dials.find_hot_pixels.debug.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0 and len(reflections) == 0:
    parser.print_help()
    exit(0)

  if len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
  if len(reflections) == 0:
    raise Sorry("No reflection lists found in input")
  if len(reflections) > 1:
    raise Sorry("Multiple reflections lists provided in input")

  assert(len(reflections) == 1)
  reflections = reflections[0]

  mask = hot_pixel_mask(imagesets[0], reflections)
  pickle.dump(mask, open(params.output.mask, 'w'), pickle.HIGHEST_PROTOCOL)

  print 'Wrote hot pixel mask to %s' % params.output.mask
  return

def hot_pixel_mask(imageset, reflections):
  depth = imageset.get_array_range()[1] - imageset.get_array_range()[0]
  xylist = filter_reflections(reflections, depth)

  from dials.array_family import flex

  mask = flex.bool(flex.grid(reversed(imageset.get_image_size())), True)

  for x, y in xylist:
    mask[y, x] = False

  print 'Found %d hot pixels' % len(xylist)

  return (mask,)

def filter_reflections(reflections, depth):
  xylist = []

  for i in range(len(reflections)):
    refl = reflections[i]
    s = refl['shoebox']
    if s.zsize() == depth:
      xylist.append((s.xoffset(), s.yoffset()))

  return xylist

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
