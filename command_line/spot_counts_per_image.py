from dials.util.options import OptionParser
from dials.util.options import flatten_reflections, flatten_datablocks
from dials.algorithms.peak_finding import per_image_analysis

import iotbx.phil
phil_scope = iotbx.phil.parse("""\
plot=None
  .type = path
individual_plots=False
  .type = bool
id = None
  .type = int(value_min=0)
""")

def run(args):
  parser = OptionParser(
    read_reflections=True,
    read_datablocks=True,
    phil=phil_scope,
    check_format=False)

  params, options = parser.parse_args(show_diff_phil=False)
  reflections = flatten_reflections(params.input.reflections)
  datablocks = flatten_datablocks(params.input.datablock)

  assert len(reflections) == 1
  assert len(datablocks) == 1

  reflections = reflections[0]
  imageset = datablocks[0].extract_imagesets()[0]

  if params.id is not None:
    reflections = reflections.select(reflections['id'] == params.id)

  stats = per_image_analysis.stats_imageset(
    imageset, reflections, plot=params.individual_plots)
  per_image_analysis.print_table(stats)
  if params.plot is not None:
    per_image_analysis.plot_stats(stats, filename=params.plot)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
