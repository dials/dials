# python imports
from __future__ import absolute_import, division
import os
from libtbx.test_utils import open_tmp_directory

def test_command_line(dials_regression):
  data_dir = os.path.join(dials_regression, 'refinement_test_data',
                          'multi_narrow_wedges')
  import glob
  experiments = glob.glob(
    os.path.join(data_dir, 'data/sweep_*/experiments.json'))

  # work in a temporary directory
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="tst_cluster_unit_cell")
  os.chdir(tmp_dir)

  import dials.util.procrunner
  result = dials.util.procrunner.run_process(
    command=['dials.cluster_unit_cell', 'plot.show=False'] + experiments,
    print_stdout=False,
    )
  assert not result['exitcode']

  #print result
  assert os.path.exists('cluster_unit_cell.png')

  from dxtbx.model.experiment_list import ExperimentList, ExperimentListFactory
  from dials.command_line import cluster_unit_cell

  experiments = ExperimentList([
    ExperimentListFactory.from_json_file(expt, check_format=False)[0]
    for expt in experiments])
  params = cluster_unit_cell.phil_scope.extract()
  params.plot.show = False
  params.plot.name = None
  clusters = cluster_unit_cell.do_cluster_analysis(experiments, params)
  assert len(clusters) == 1
  cluster = clusters[0]
  assert len(cluster.members) == 40
  from libtbx.test_utils import approx_equal
  assert approx_equal(
    cluster.medians,
    [90.9430182020995, 90.9430182020995, 90.9430182020995,
     109.47122063449069, 109.47122063449069, 109.47122063449069])
  assert approx_equal(
    cluster.stdevs,
    [0.09509739126548639, 0.09509739126548526, 0.0950973912654865, 0, 0, 0])
