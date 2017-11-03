from __future__ import absolute_import, division

import iotbx.phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments

help_message = '''
Augment spot list with additional information - for example number of pixels
in peak region etc.'''

phil_scope = iotbx.phil.parse('''
output {
  reflections = stronger.pickle
  .type = path
}
''')

def augment_reflections(reflections, datablock=None):
  '''Add extra columns of derived data.'''

  from dials.algorithms.shoebox import MaskCode
  good = (MaskCode.Foreground | MaskCode.Valid)

  x0, x1, y0, y1, z0, z2 = reflections['bbox'].parts()

  dx = x1 - x0
  dy = y1 - y0

  # compute signal pixels in each shoebox as an array
  n_signal = reflections['shoebox'].count_mask_values(good)

  reflections['dx'] = dx
  reflections['dy'] = dy
  reflections['n_signal'] = n_signal

  return reflections

def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  from dials.util import log
  usage = "%s [options] [datablock.json] strong.pickle" % \
    libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)

  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

  if len(reflections) != 1:
    raise Sorry("Exactly one reflection file needed")
  if len(datablocks) > 1:
    raise Sorry("0, 1 datablocks required")

  datablock = None
  if len(datablocks) == 1:
    datablock = datablocks[0]

  stronger = augment_reflections(reflections[0], datablock=datablock)
  stronger.as_pickle(params.output.reflections)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
