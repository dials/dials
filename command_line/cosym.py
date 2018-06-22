from __future__ import absolute_import, division, print_function
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


import logging
logger = logging.getLogger('dials.command_line.cosym')

import os
from libtbx.utils import Sorry
import iotbx.phil
from cctbx import crystal, miller
from cctbx import sgtbx
from iotbx.reflection_file_reader import any_reflection_file
from dials.array_family import flex
from dials.util.options import flatten_experiments, flatten_reflections
from dials.algorithms.symmetry.cosym import analyse_datasets

phil_scope = iotbx.phil.parse('''\
d_min = Auto
  .type = float(value_min=0)

min_i_mean_over_sigma_mean = None
  .type = float(value_min=0)

batch = None
  .type = ints(value_min=0, size=2)

normalisation = kernel
  .type = choice

mode = *full ambiguity
  .type = choice

space_group = None
  .type = space_group

unit_cell_clustering {
  threshold = 5000
    .type = float(value_min=0)
    .help = 'Threshold value for the clustering'
  log = False
    .type = bool
    .help = 'Display the dendrogram with a log scale'
}

include scope dials.algorithms.symmetry.cosym.phil_scope

seed = 230
  .type = int(value_min=0)

output {
  suffix = "_reindexed"
    .type = str
  log = cosym.log
    .type = str
  debug_log = cosym.debug.log
    .type = str
  experiments = "reindexed_experiments.json"
    .type = path
  reflections = "reindexed_reflections.pickle"
    .type = path
}

verbosity = 1
  .type = int(value_min=0)
  .help = "The verbosity level"
''', process_includes=True)

def run(args):
  import libtbx
  from libtbx import easy_pickle
  from dials.util import log
  from dials.util.options import OptionParser

  parser = OptionParser(
    #usage=usage,
    phil=phil_scope,
    read_reflections=True,
    read_datablocks=False,
    read_experiments=True,
    check_format=False,
    #epilog=help_message
  )

  params, options, args = parser.parse_args(
    args=args, show_diff_phil=False, return_unhandled=True)

  # Configure the logging
  log.config(
    params.verbosity,
    info=params.output.log,
    debug=params.output.debug_log)

  from dials.util.version import dials_version
  logger.info(dials_version())

  # Log the diff phil
  diff_phil = parser.diff_phil.as_str()
  if diff_phil is not '':
    logger.info('The following parameters have been modified:\n')
    logger.info(diff_phil)

  if params.seed is not None:
    import random
    flex.set_random_seed(params.seed)
    random.seed(params.seed)

  if params.save_plot and not params.animate:
    import matplotlib
    # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
    matplotlib.use('Agg') # use a non-interactive backend

  datasets_input = []

  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) or len(reflections):
    if len(reflections) == 1:
      reflections_input = reflections[0]
      reflections = []
      for i in range(len(experiments)):
        reflections.append(reflections_input.select(reflections_input['id'] == i))

    if len(experiments) > len(reflections):
      flattened_reflections = []
      for refl in reflections:
        for i in range(0, flex.max(refl['id'])+1):
          sel = refl['id'] == i
          flattened_reflections.append(refl.select(sel))
      reflections = flattened_reflections

    assert len(experiments) == len(reflections)

    i_refl = 0
    for i_expt in enumerate(experiments):
      refl = reflections[i_refl]

    for expt, refl in zip(experiments, reflections):
      crystal_symmetry = crystal.symmetry(
        unit_cell=expt.crystal.get_unit_cell(),
        space_group=expt.crystal.get_space_group())
      if 0 and 'intensity.prf.value' in refl:
        sel = refl.get_flags(refl.flags.integrated_prf)
        assert sel.count(True) > 0
        refl = refl.select(sel)
        data = refl['intensity.prf.value']
        variances = refl['intensity.prf.variance']
      else:
        assert 'intensity.sum.value' in refl
        sel = refl.get_flags(refl.flags.integrated_sum)
        assert sel.count(True) > 0
        refl = refl.select(sel)
        data = refl['intensity.sum.value']
        variances = refl['intensity.sum.variance']
      # FIXME probably need to do some filtering of intensities similar to that
      # done in export_mtz
      miller_indices = refl['miller_index']
      assert variances.all_gt(0)
      sigmas = flex.sqrt(variances)

      miller_set = miller.set(crystal_symmetry, miller_indices, anomalous_flag=False)
      intensities = miller.array(miller_set, data=data, sigmas=sigmas)
      intensities.set_observation_type_xray_intensity()
      intensities.set_info(miller.array_info(
        source='DIALS',
        source_type='pickle'
      ))
      datasets_input.append(intensities)

  files = args

  for file_name in files:

    try:
      data = easy_pickle.load(file_name)
      intensities = data['observations'][0]
      intensities.set_info(miller.array_info(
        source=file_name,
        source_type='pickle'
      ))
      intensities = intensities.customized_copy(
        anomalous_flag=False).set_info(intensities.info())
      batches = None
    except Exception:
      reader = any_reflection_file(file_name)
      assert reader.file_type() == 'ccp4_mtz'

      as_miller_arrays = reader.as_miller_arrays(merge_equivalents=False)
      intensities = [ma for ma in as_miller_arrays
                     if ma.info().labels == ['I', 'SIGI']][0]
      batches = [ma for ma in as_miller_arrays
                 if ma.info().labels == ['BATCH']]
      if len(batches):
        batches = batches[0]
      else:
        batches = None
      mtz_object = reader.file_content()
      intensities = intensities.customized_copy(
        anomalous_flag=False,
        indices=mtz_object.extract_original_index_miller_indices()).set_info(
          intensities.info())

    intensities.set_observation_type_xray_intensity()
    datasets_input.append(intensities)

  if len(datasets_input) == 0:
    raise Sorry('No valid reflection files provided on command line')

  datasets = []

  # per-dataset change of basis operator to ensure all consistent
  change_of_basis_ops = []

  for intensities in datasets_input:

    if params.batch is not None:
      assert batches is not None
      bmin, bmax = params.batch
      assert bmax >= bmin
      sel = (batches.data() >= bmin) & (batches.data() <= bmax)
      assert sel.count(True) > 0
      intensities = intensities.select(sel)

    if params.min_i_mean_over_sigma_mean is not None and (
         params.d_min is libtbx.Auto or params.d_min is not None):
      from xia2.Modules import Resolutionizer
      rparams = Resolutionizer.phil_defaults.extract().resolutionizer
      rparams.nbins = 20
      resolutionizer = Resolutionizer.resolutionizer(intensities, rparams, None)
      i_mean_over_sigma_mean = 4
      d_min = resolutionizer.resolution_i_mean_over_sigma_mean(i_mean_over_sigma_mean)
      if params.d_min is libtbx.Auto:
        intensities = intensities.resolution_filter(d_min=d_min).set_info(
          intensities.info())
        if params.verbose:
          logger.info('Selecting reflections with d > %.2f' %d_min)
      elif d_min > params.d_min:
        logger.info(
          'Rejecting dataset %s as d_min too low (%.2f)' %(file_name, d_min))
        continue
      else:
        logger.info('Estimated d_min for %s: %.2f' %(file_name, d_min))
    elif params.d_min not in (None, libtbx.Auto):
      intensities = intensities.resolution_filter(d_min=params.d_min).set_info(
        intensities.info())

    if params.normalisation == 'kernel':
      from mmtbx.scaling import absolute_scaling
      normalisation = absolute_scaling.kernel_normalisation(intensities, auto_kernel=True)
      intensities = normalisation.normalised_miller.deep_copy()

    cb_op_to_primitive = intensities.change_of_basis_op_to_primitive_setting()
    change_of_basis_ops.append(cb_op_to_primitive)
    intensities = intensities.change_basis(cb_op_to_primitive)
    if params.mode == 'full' or params.space_group is not None:
      if params.space_group is not None:
        space_group_info = params.space_group.primitive_setting()
        if not space_group_info.group().is_compatible_unit_cell(intensities.unit_cell()):
          logger.info(
            'Skipping data set - incompatible space group and unit cell: %s, %s' %(
              space_group_info, intensities.unit_cell()))
          continue
      else:
        space_group_info = sgtbx.space_group_info('P1')
      intensities = intensities.customized_copy(space_group_info=space_group_info)

    datasets.append(intensities)

  crystal_symmetries = [d.crystal_symmetry().niggli_cell() for d in datasets]
  lattice_ids = range(len(datasets))
  from xfel.clustering.cluster import Cluster
  from xfel.clustering.cluster_groups import unit_cell_info
  ucs = Cluster.from_crystal_symmetries(crystal_symmetries, lattice_ids=lattice_ids)
  threshold = 1000
  if params.save_plot:
    from matplotlib import pyplot as plt
    fig = plt.figure("Andrews-Bernstein distance dendogram", figsize=(12, 8))
    ax = plt.gca()
  else:
    ax = None
  clusters, _ = ucs.ab_cluster(
    params.unit_cell_clustering.threshold,
    log=params.unit_cell_clustering.log,
    write_file_lists=False,
    schnell=False,
    doplot=params.save_plot,
    ax=ax
  )
  if params.save_plot:
    plt.tight_layout()
    plt.savefig('%scluster_unit_cell.png' % params.plot_prefix)
    plt.close(fig)
  logger.info(unit_cell_info(clusters))
  largest_cluster = None
  largest_cluster_lattice_ids = None
  for cluster in clusters:
    cluster_lattice_ids = [m.lattice_id for m in cluster.members]
    if largest_cluster_lattice_ids is None:
      largest_cluster_lattice_ids = cluster_lattice_ids
    elif len(cluster_lattice_ids) > len(largest_cluster_lattice_ids):
      largest_cluster_lattice_ids = cluster_lattice_ids

  dataset_selection = largest_cluster_lattice_ids
  if len(dataset_selection) < len(datasets):
    logger.info(
      'Selecting subset of data for cosym analysis: %s' %str(dataset_selection))
    datasets = [datasets[i] for i in dataset_selection]

  for i, dataset in enumerate(datasets):
    metric_subgroups = sgtbx.lattice_symmetry.metric_subgroups(dataset, max_delta=5)
    subgroup = metric_subgroups.result_groups[0]
    cb_op_inp_best = subgroup['cb_op_inp_best']
    datasets[i] = dataset.change_basis(cb_op_inp_best)
    change_of_basis_ops[i] = cb_op_inp_best * change_of_basis_ops[i]

  cb_op_ref_min = datasets[0].change_of_basis_op_to_niggli_cell()
  for i, dataset in enumerate(datasets):
    if params.space_group is None:
      datasets[i] = dataset.change_basis(cb_op_ref_min).customized_copy(
        space_group_info=sgtbx.space_group_info('P1'))
    else:
      datasets[i] = dataset.change_basis(cb_op_ref_min)
      datasets[i] = datasets[i].customized_copy(
        crystal_symmetry=crystal.symmetry(
          unit_cell=datasets[i].unit_cell(),
          space_group_info=params.space_group.primitive_setting(),
          assert_is_compatible_unit_cell=False))
    datasets[i] = datasets[i].merge_equivalents().array()
    change_of_basis_ops[i] = cb_op_ref_min * change_of_basis_ops[i]

  result = analyse_datasets(datasets, params)

  space_groups = {}
  reindexing_ops = {}
  for dataset_id in result.reindexing_ops.iterkeys():
    if 0 in result.reindexing_ops[dataset_id]:
      cb_op = result.reindexing_ops[dataset_id][0]
      reindexing_ops.setdefault(cb_op, [])
      reindexing_ops[cb_op].append(dataset_id)
    if dataset_id in result.space_groups:
      space_groups.setdefault(result.space_groups[dataset_id], [])
      space_groups[result.space_groups[dataset_id]].append(dataset_id)

  logger.info('Space groups:')
  for sg, datasets in space_groups.iteritems():
    logger.info(str(sg.info().reference_setting()))
    logger.info(datasets)

  logger.info('Reindexing operators:')
  for cb_op, datasets in reindexing_ops.iteritems():
    logger.info(cb_op)
    logger.info(datasets)

  if (len(experiments) and len(reflections) and
      params.output.reflections is not None and
      params.output.experiments is not None):
    import copy
    from dxtbx.model import ExperimentList
    from dxtbx.serialize import dump
    reindexed_experiments = ExperimentList()
    reindexed_reflections = flex.reflection_table()
    expt_id = 0
    for cb_op, dataset_ids in reindexing_ops.iteritems():
      cb_op = sgtbx.change_of_basis_op(cb_op)
      for dataset_id in dataset_ids:
        expt = experiments[dataset_selection[dataset_id]]
        refl = reflections[dataset_selection[dataset_id]]
        reindexed_expt = copy.deepcopy(expt)
        refl_reindexed = copy.deepcopy(refl)
        cb_op_this = cb_op * change_of_basis_ops[dataset_id]
        reindexed_expt.crystal = reindexed_expt.crystal.change_basis(cb_op_this)
        refl_reindexed['miller_index'] = cb_op_this.apply(
          refl_reindexed['miller_index'])
        reindexed_experiments.append(reindexed_expt)
        refl_reindexed['id'] = flex.int(refl_reindexed.size(), expt_id)
        reindexed_reflections.extend(refl_reindexed)
        expt_id += 1

    logger.info('Saving reindexed experiments to %s' % params.output.experiments)
    dump.experiment_list(reindexed_experiments, params.output.experiments)
    logger.info('Saving reindexed reflections to %s' % params.output.reflections)
    reindexed_reflections.as_pickle(params.output.reflections)

  elif params.output.suffix is not None:
    for cb_op, dataset_ids in reindexing_ops.iteritems():
      cb_op = sgtbx.change_of_basis_op(cb_op)
      for dataset_id in dataset_ids:
        file_name = files[dataset_selection[dataset_id]]
        basename = os.path.basename(file_name)
        out_name = os.path.splitext(
          basename)[0] + params.output.suffix + '_' + str(dataset_selection[dataset_id]) + ".mtz"
        reader = any_reflection_file(file_name)
        assert reader.file_type() == 'ccp4_mtz'
        mtz_object = reader.file_content()
        cb_op_this = cb_op * change_of_basis_ops[dataset_id]
        if not cb_op_this.is_identity_op():
          logger.info('reindexing %s (%s)' %(file_name, cb_op_this.as_xyz()))
          mtz_object.change_basis_in_place(cb_op_this)
        mtz_object.write(out_name)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
