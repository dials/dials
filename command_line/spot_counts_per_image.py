from __future__ import division

from libtbx import phil
from dials.util.options import OptionParser
from dials.util.options import flatten_reflections
from dials.util.options import flatten_datablocks
from dials.array_family import flex

master_phil_scope = phil.parse("""\
output {
  plot = None
    .type = path
}
""")

def run(args):
  parser = OptionParser(
    phil=master_phil_scope,
    read_reflections=True,
    read_datablocks=True,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)
  if len(datablocks) == 0:
    print "No Datablock could be constructed"
    return
  elif len(datablocks) > 1:
    raise RuntimeError("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
  assert(len(reflections) == 1)
  reflections = reflections[0]

  image_number = reflections['xyzobs.px.value'].parts()[2]
  image_number = flex.floor(image_number)

  hist = flex.histogram(
    image_number, data_min=0, data_max=len(imagesets[0]),
    n_slots=len(imagesets[0]))

  print "Per-image analysis:"
  print_table(hist)

  if params.output.plot is not None:
    plot(hist, params.output.plot)

  return

def table(histogram):
  rows = [("image", "#spots")]
  for i_image, spot_count in zip(
      histogram.slot_centers(), histogram.slots()):
    rows.append((str(int(i_image+0.5)), str(spot_count)))
  return rows

def print_table(histogram, out=None):
  if out is None:
    import sys
    out = sys.stdout
  from libtbx import table_utils

  rows = table(histogram)
  print >> out, table_utils.format(
    rows, has_header=True, prefix="|", postfix="|")


def plot(histogram, file_name):
  try:
    import matplotlib
    matplotlib.use('Agg') # use a non-interactive backend
    # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
    from matplotlib import pyplot
  except ImportError:
    raise Sorry("matplotlib must be installed to generate a plot.")

  n_spots = histogram.slots()
  #resolution = columns.get('resolution')
  i_image = histogram.slot_centers() + 0.5
  fig = pyplot.figure()
  ax1 = fig.add_subplot(111)
  sc1 = ax1.scatter(i_image, n_spots, s=20, color='blue', marker='o', alpha=0.5)
  ax1.set_xlabel('Image #')
  ax1.set_ylabel('# spots')
  ax1.set_xlim((0.0, len(n_spots)))
  ax1.set_ylim(bottom=-0.2)
  #ax2 = ax1.twinx()
  #sc2 = ax2.scatter(i_image, resolution, s=20, color='red', marker='^', alpha=0.5)
  #ax2.set_ylabel(u'resolution (\u00c5)')
  #ax2.set_xlim((0, len(n_spots)))
  #ax2.invert_yaxis()

  # Use mode="fixed" as mode="expand" causes floating point error on some
  # versions of matplotlib.
  # See https://github.com/matplotlib/matplotlib/pull/1864
  lgd = pyplot.legend(
    (sc1, ), ('#spots',), ncol=1,
    loc='upper center',
    mode="fixed", borderaxespad=0.,
    bbox_to_anchor=(0.0,-0.22, 1., .102))
  pyplot.savefig(file_name, dpi=600, bbox_extra_artists=(lgd,),
                 bbox_inches='tight')



if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
