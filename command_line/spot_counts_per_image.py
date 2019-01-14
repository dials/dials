from __future__ import absolute_import, division, print_function
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

import sys

from dials.util.options import OptionParser
from dials.util.options \
     import flatten_reflections, flatten_datablocks, flatten_experiments
from dials.algorithms.spot_finding import per_image_analysis

import iotbx.phil

help_message = '''

Reports the number of strong spots and computes an estimate of the resolution
limit for each image, given the results of dials.find_spots. Optionally
generates a plot of the per-image statistics (plot=image.png).

Examples::

  dials.spot_counts_per_image datablock.json strong.pickle

  dials.spot_counts_per_image datablock.json strong.pickle plot=per_image.png

'''

phil_scope = iotbx.phil.parse("""\
resolution_analysis = True
  .type = bool
plot = None
  .type = path
json = None
  .type = path
split_json = False
  .type = bool
joint_json = True
  .type = bool
individual_plots = False
  .type = bool
id = None
  .type = int(value_min=0)
""")

def run(args):
  usage = "dials.spot_counts_per_image [options] datablock.json strong.pickle"

  parser = OptionParser(
    usage=usage,
    read_reflections=True,
    read_datablocks=True,
    read_experiments=True,
    phil=phil_scope,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)
  reflections = flatten_reflections(params.input.reflections)
  datablocks = flatten_datablocks(params.input.datablock)
  experiments = flatten_experiments(params.input.experiments)

  if not any([reflections, experiments, datablocks]):
    parser.print_help()
    return

  if len(reflections) != 1:
    sys.exit('exactly 1 reflection table must be specified')
  if len(datablocks) != 1:
    if experiments:
      if len(experiments.imagesets()) != 1:
        sys.exit('exactly 1 datablock must be specified')
      imageset = experiments.imagesets()[0]
    else:
      sys.exit('exactly 1 datablock must be specified')
  else:
    imageset = datablocks[0].extract_imagesets()[0]

  reflections = reflections[0]

  if params.id is not None:
    reflections = reflections.select(reflections['id'] == params.id)

  stats = per_image_analysis.stats_imageset(
    imageset, reflections, resolution_analysis=params.resolution_analysis,
    plot=params.individual_plots)
  per_image_analysis.print_table(stats)

  from libtbx import table_utils
  overall_stats = per_image_analysis.stats_single_image(
    imageset, reflections, resolution_analysis=params.resolution_analysis)
  rows = [
    ("Overall statistics", ""),
    ("#spots", "%i" % overall_stats.n_spots_total),
    ("#spots_no_ice", "%i" % overall_stats.n_spots_no_ice),
    ("d_min", "%.2f" % overall_stats.estimated_d_min),
    ("d_min (distl method 1)", "%.2f (%.2f)" % (
      overall_stats.d_min_distl_method_1, overall_stats.noisiness_method_1)),
    ("d_min (distl method 2)", "%.2f (%.2f)" % (
      overall_stats.d_min_distl_method_1, overall_stats.noisiness_method_1)),
    ]
  print(table_utils.format(rows, has_header=True, prefix="| ", postfix=" |"))

  if params.json:
    import json
    if params.split_json:
      for k in stats.__dict__:
        start, end = params.json.split('.')
        with open('%s_%s.%s' % (start, k, end), 'wb') as fp:
          json.dump(stats.__dict__[k], fp)
    if params.joint_json:
      with open(params.json, 'wb') as fp:
        json.dump(stats.__dict__, fp)
  if params.plot:
    per_image_analysis.plot_stats(stats, filename=params.plot)

if __name__ == '__main__':
  run(sys.argv[1:])
