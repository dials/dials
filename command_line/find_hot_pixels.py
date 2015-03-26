from __future__ import division

from libtbx.phil import command_line
import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments

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
Find hot pixels from custom dials.find_spots run, i.e.

  dials.find_hot_pixels datablock.json strong.pickle

Will return a hot_pixel.pickle mask for use in integration.
'''

def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util import log
  from logging import info
  import cPickle as pickle
  usage = "%s [options] datablock.json strong.pickle" % \
    libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  log.config(
    params.verbosity, info='dials.find_hot_pixels.log',
      debug='dials.find_hot_pixels.debug.log')

  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

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
