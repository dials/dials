# LIBTBX_SET_DISPATCHER_NAME dev.dials.filter_good_intensities

from __future__ import division
from dials.util.export_mtz import export_mtz

def filter_good_reflections(integrated_data):
  assert(min(integrated_data['id']) == max(integrated_data['id']) == 0)

  selection = integrated_data['intensity.sum.variance'] <= 0
  if selection.count(True) > 0:
    integrated_data.del_selected(selection)
    print 'Removing %d reflections with negative variance' % \
          selection.count(True)

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      print 'Removing %d profile reflections with negative variance' % \
            selection.count(True)

  if 'partiality' in integrated_data:
    selection = integrated_data['partiality'] < 0.99
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      print 'Removing %d incomplete reflections' % \
        selection.count(True)

  return integrated_data

def run(args):

  import libtbx.load_env
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  from libtbx.utils import Sorry
  from libtbx.phil import parse
  phil_scope = parse('''
    hklout = hklout.pickle
      .type = str
      .help = "The output pickle file"
  ''')

  usage = '%s integrated.pickle [hklout=hklout.pickle]' % (
              libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage = usage,
    read_reflections=True,
    check_format=False,
    phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  reflections = flatten_reflections(params.input.reflections)
  if len(reflections) != 1:
    raise Sorry('exactly 1 reflection table must be specified')

  integrated_data = reflections[0]
  filter_good_reflections(integrated_data).as_pickle(params.hklout)




if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
