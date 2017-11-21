# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1

from __future__ import absolute_import, division
import iotbx.phil
from cctbx import sgtbx

help_message = '''
'''

phil_scope= iotbx.phil.parse('''
threshold = 5000
  .type = float(value_min=0)
  .help = 'Threshold value for the clustering'
plot {
  show = True
    .type = bool
  name = 'cluster_unit_cell.png'
    .type = path
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}
''')


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  import libtbx.load_env

  usage = "%s [options] datablock.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)

  if len(experiments) == 0:
    parser.print_help()
    exit(0)

  do_cluster_analysis(experiments, params)

def do_cluster_analysis(experiments, params):

  from cctbx import crystal
  from xfel.clustering.cluster import Cluster
  from xfel.clustering.cluster_groups import unit_cell_info

  crystal_symmetries = [
    crystal.symmetry(unit_cell=expt.crystal.get_unit_cell(),
                     space_group=expt.crystal.get_space_group())
    for expt in experiments]
  ucs = Cluster.from_crystal_symmetries(crystal_symmetries)

  if params.plot.show or params.plot.name is not None:
    if not params.plot.show:
      import matplotlib
      # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
      matplotlib.use('Agg') # use a non-interactive backend
    import matplotlib.pyplot as plt
    plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
    ax = plt.gca()
    clusters, cluster_axes = ucs.ab_cluster(
      params.threshold,
      log=params.plot.log,
      ax=ax,
      write_file_lists=False,
      #schnell=_args.schnell,
      doplot=True)
    print unit_cell_info(clusters)
    plt.tight_layout()
    if params.plot.name is not None:
      plt.savefig(params.plot.name)
    if params.plot.show:
      plt.show()

  else:
    clusters, cluster_axes = ucs.ab_cluster(
      params.threshold,
      log=params.plot.log,
      write_file_lists=False,
      #schnell=_args.schnell,
      doplot=False)
    print unit_cell_info(clusters)

  return clusters


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
