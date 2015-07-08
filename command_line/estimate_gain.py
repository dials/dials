from __future__ import division

import iotbx.phil
from scitbx.array_family import flex
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks

help_message = '''

'''

phil_scope = iotbx.phil.parse("""\
kernel_size = 10,10
  .type = ints(size=2, value_min=1)
""", process_includes=True)

def estimate_gain(imageset, kernel_size=(10,10)):
  detector = imageset.get_detector()

  from dials.algorithms.image.threshold import KabschDebug

  raw_data = imageset[0]
  if len(detector) == 1:
    raw_data = [raw_data]

  gain_value = 1
  gain_map = [flex.double(raw_data[i].accessor(), gain_value)
              for i in range(len(detector))]

  mask = imageset.get_mask(0)

  min_local = 0

  # dummy values, shouldn't affect results
  nsigma_b = 6
  nsigma_s = 3
  global_threshold = 0

  kabsch_debug_list = []
  for i_panel in range(len(detector)):
    kabsch_debug_list.append(
      KabschDebug(
        raw_data[i_panel].as_double(), mask[i_panel], gain_map[i_panel],
        kernel_size, nsigma_b, nsigma_s, global_threshold, min_local))

  dispersion = flex.double()
  for kabsch in kabsch_debug_list:
    dispersion.extend(kabsch.coefficient_of_variation().as_1d())

  sorted_dispersion = flex.sorted(dispersion)
  from libtbx.math_utils import nearest_integer as nint

  q1 = sorted_dispersion[nint(len(sorted_dispersion)/4)]
  q2 = sorted_dispersion[nint(len(sorted_dispersion)/2)]
  q3 = sorted_dispersion[nint(len(sorted_dispersion)*3/4)]
  iqr = q3-q1

  print "q1, q2, q3: %.2f, %.2f, %.2f" %(q1, q2, q3)

  inlier_sel = (sorted_dispersion > (q1 - 1.5*iqr)) & (sorted_dispersion < (q3 + 1.5*iqr))
  sorted_dispersion = sorted_dispersion.select(inlier_sel)
  gain = sorted_dispersion[nint(len(sorted_dispersion)/2)]
  print "Estimated gain: %.2f" %gain

  if 0:
    sel = flex.random_selection(population_size=len(sorted_dispersion), sample_size=10000)
    sorted_dispersion = sorted_dispersion.select(sel)

    from matplotlib import pyplot
    pyplot.scatter(range(len(sorted_dispersion)), sorted_dispersion)
    pyplot.ylim(0, 10)
    pyplot.show()

  return gain


def run(args):
  import libtbx.load_env
  from libtbx.utils import Sorry
  usage = "%s [options] datablock.json" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    check_format=True,
    read_datablocks_from_images=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)

  ## Configure the logging
  #log.config(
    #params.verbosity, info='dials.estimate_gain.log', debug='dials.estimate_gain.debug.log')

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    print 'The following parameters have been modified:\n'
    print diff_phil

  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    return
  elif len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = []
    for datablock in datablocks:
      imagesets.extend(datablock.extract_imagesets())

  assert len(imagesets) == 1
  imageset = imagesets[0]
  estimate_gain(imageset, params.kernel_size)

  return


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
