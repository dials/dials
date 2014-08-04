from __future__ import division
import os
import math
import glob

import libtbx
from libtbx import easy_mp, easy_pickle
from libtbx import group_args
from libtbx import table_utils
from dials.array_family import flex
from cctbx import crystal

def run_once(directory):
  from dxtbx.serialize import load
  sweep_dir = os.path.basename(directory)
  print sweep_dir

  datablock_name = os.path.join(directory, "datablock.json")
  if not os.path.exists(datablock_name):
    # this is what xia2 calls it:
    datablock_name = os.path.join(directory, "datablock_import.json")
  strong_spots_name = os.path.join(directory, "strong.pickle")
  experiments_name = os.path.join(directory, "experiments.json")
  indexed_spots_name = os.path.join(directory, "indexed.pickle")
  unindexed_spots_name = os.path.join(directory, "unindexed.pickle")
  if not (os.path.exists(datablock_name) and os.path.exists(strong_spots_name)):
    return

  datablock = load.datablock(datablock_name)
  assert len(datablock) == 1
  if len(datablock[0].extract_sweeps()) == 0:
    print "Skipping %s" %directory
    return
  sweep = datablock[0].extract_sweeps()[0]
  template = sweep.get_template()

  strong_spots = easy_pickle.load(strong_spots_name)
  n_strong_spots = len(strong_spots)
  if os.path.exists(experiments_name):
    experiments = load.experiment_list(experiments_name)
    n_indexed_lattices = len(experiments)
  else:
    experiments = None
    n_indexed_lattices = 0

  g = glob.glob(os.path.join(directory, "xds*", "run_2", "INTEGRATE.HKL"))
  n_integrated_lattices = len(g)

  if os.path.exists(indexed_spots_name):
    indexed_spots = easy_pickle.load(indexed_spots_name)
  else:
    indexed_spots = None
    g = glob.glob(os.path.join(directory, "indexed_*.pickle"))
    if len(g):
      for path in g:
        if indexed_spots is None:
          indexed_spots = easy_pickle.load(path)
        else:
          indexed_spots.extend(easy_pickle.load(path))

  if os.path.exists(unindexed_spots_name):
    unindexed_spots = easy_pickle.load(unindexed_spots_name)
    n_unindexed_spots = len(unindexed_spots)
  else:
    n_unindexed_spots = 0

  # calculate estimated d_min for sweep based on 95th percentile
  from dials.algorithms.indexing import indexer2
  detector = sweep.get_detector()
  scan = sweep.get_scan()
  beam = sweep.get_beam()
  goniometer = sweep.get_goniometer()
  if len(strong_spots) == 0:
    d_strong_spots_99th_percentile = 0
    d_strong_spots_95th_percentile = 0
    d_strong_spots_50th_percentile = 0
    n_strong_spots_dmin_4 = 0
  else:
    spots_mm = indexer2.indexer_base.map_spots_pixel_to_mm_rad(
      strong_spots, detector, scan)
    indexer2.indexer_base.map_centroids_to_reciprocal_space(
      spots_mm, detector, beam, goniometer)
    d_spacings = 1/spots_mm['rlp'].norms()
    perm = flex.sort_permutation(d_spacings, reverse=True)
    d_spacings_sorted = d_spacings.select(perm)
    percentile_99th = int(math.floor(0.99 * len(d_spacings)))
    percentile_95th = int(math.floor(0.95 * len(d_spacings)))
    percentile_50th = int(math.floor(0.5 * len(d_spacings)))
    d_strong_spots_99th_percentile = d_spacings_sorted[percentile_99th]
    d_strong_spots_95th_percentile = d_spacings_sorted[percentile_95th]
    d_strong_spots_50th_percentile = d_spacings_sorted[percentile_50th]
    n_strong_spots_dmin_4 = (d_spacings >= 4).count(True)

  cell_params = flex.sym_mat3_double()
  n_indexed = flex.double()
  d_min_indexed = flex.double()
  rmsds = flex.vec3_double()
  sweep_dir_cryst = flex.std_string()
  if experiments is not None:
    for i, experiment in enumerate(experiments):
      sweep_dir_cryst.append(sweep_dir)
      crystal_model = experiment.crystal
      unit_cell = crystal_model.get_unit_cell()
      space_group = crystal_model.get_space_group()
      crystal_symmetry = crystal.symmetry(unit_cell=unit_cell,
                                          space_group=space_group)
      cb_op_reference_setting =  crystal_symmetry.change_of_basis_op_to_reference_setting()
      crystal_symmetry_reference_setting = crystal_symmetry.change_basis(
        cb_op_reference_setting)
      cell_params.append(crystal_symmetry_reference_setting.unit_cell().parameters())
      spots_mm = indexed_spots.select(indexed_spots['id'] == i)
      n_indexed.append(len(spots_mm))
      if len(spots_mm) == 0:
        d_min_indexed.append(0)
      else:
        indexer2.indexer_base.map_centroids_to_reciprocal_space(
          spots_mm, detector, beam, goniometer)
        d_spacings = 1/spots_mm['rlp'].norms()
        perm = flex.sort_permutation(d_spacings, reverse=True)
        d_min_indexed.append(d_spacings[perm[-1]])
      try:
        rmsds.append(get_rmsds_obs_pred(spots_mm, experiment))
      except Exception, e:
        print e
        rmsds.append((-1,-1,-1))
        continue

  return group_args(
    sweep_dir=sweep_dir,
    template=template,
    n_strong_spots=n_strong_spots,
    n_strong_spots_dmin_4=n_strong_spots_dmin_4,
    n_unindexed_spots=n_unindexed_spots,
    n_indexed_lattices=n_indexed_lattices,
    n_integrated_lattices=n_integrated_lattices,
    d_strong_spots_50th_percentile=d_strong_spots_50th_percentile,
    d_strong_spots_95th_percentile=d_strong_spots_95th_percentile,
    d_strong_spots_99th_percentile=d_strong_spots_99th_percentile,
    cell_params=cell_params,
    n_indexed=n_indexed,
    d_min_indexed=d_min_indexed,
    rmsds=rmsds,
    sweep_dir_cryst=sweep_dir_cryst)


def get_rmsds_obs_pred(observations, experiment):
  from dials.algorithms.spot_prediction import ray_intersection
  from dials.algorithms.indexing.indexer2 import master_params
  from dials.algorithms.refinement import RefinerFactory
  from dxtbx.model.experiment.experiment_list import ExperimentList
  master_params.refinement.reflections.close_to_spindle_cutoff = 0.001
  from dials.model.data import ReflectionList
  ref_list = ReflectionList.from_table(observations)
  ref_list = ray_intersection(experiment.detector, ref_list)
  ref_table = ref_list.to_table()
  import copy
  reflections = copy.deepcopy(observations)
  reflections['xyzcal.mm'] = ref_table['xyzcal.mm']
  reflections['xyzcal.px'] = ref_table['xyzcal.px']

  # XXX hack to make it work for a single lattice
  reflections['id'] = flex.int(len(reflections), 0)
  refine = RefinerFactory.from_parameters_data_experiments(
    master_params, reflections, ExperimentList([experiment]), verbosity=0)
  return refine.rmsds()


def run(args):
  sweep_directories = []
  templates = []
  n_strong_spots = flex.int()
  n_strong_spots_dmin_4 = flex.int()
  d_strong_spots_99th_percentile = flex.double()
  d_strong_spots_95th_percentile = flex.double()
  d_strong_spots_50th_percentile = flex.double()
  n_unindexed_spots = flex.int()
  n_indexed_lattices = flex.int()
  n_integrated_lattices = flex.int()
  sweep_dir_cryst = flex.std_string()

  orig_dir = os.path.abspath(os.curdir)

  rmsds = flex.vec3_double()
  cell_params = flex.sym_mat3_double()
  n_indexed = flex.double()
  d_min_indexed = flex.double()
  rmsds = flex.vec3_double()

  nproc = easy_mp.get_processes(libtbx.Auto)
  #nproc = 1
  results = easy_mp.parallel_map(
    func=run_once,
    iterable=args,
    processes=nproc,
    method="multiprocessing",
    preserve_order=True,
    asynchronous=True,
    preserve_exception_message=True,
  )

  for result in results:
    if result is None: continue
    sweep_directories.append(result.sweep_dir)
    templates.append(result.template)
    n_strong_spots.append(result.n_strong_spots)
    n_strong_spots_dmin_4.append(result.n_strong_spots_dmin_4)
    n_unindexed_spots.append(result.n_unindexed_spots)
    n_indexed_lattices.append(result.n_indexed_lattices)
    n_integrated_lattices.append(result.n_integrated_lattices)
    d_strong_spots_50th_percentile.append(result.d_strong_spots_50th_percentile)
    d_strong_spots_95th_percentile.append(result.d_strong_spots_95th_percentile)
    d_strong_spots_99th_percentile.append(result.d_strong_spots_99th_percentile)
    cell_params.extend(result.cell_params)
    n_indexed.extend(result.n_indexed)
    d_min_indexed.extend(result.d_min_indexed)
    rmsds.extend(result.rmsds)
    sweep_dir_cryst.extend(result.sweep_dir_cryst)

  table_data = [('sweep_dir', 'template', '#strong_spots', '#unindexed_spots', '#lattices',
                 'd_spacing_50th_percentile', 'd_spacing_95th_percentile',
                 'd_spacing_99th_percentile',)]
  for i in range(len(sweep_directories)):
    table_data.append((sweep_directories[i],
                       templates[i],
                       str(n_strong_spots[i]),
                       str(n_unindexed_spots[i]),
                       str(n_indexed_lattices[i]),
                       str(d_strong_spots_50th_percentile[i]),
                       str(d_strong_spots_95th_percentile[i]),
                       str(d_strong_spots_99th_percentile[i]),
                       ))

  with open('results.txt', 'wb') as f:
    print >> f, table_utils.format(
      table_data, has_header=True, justify='right')

  table_data = [('sweep_dir', 'cell_a', 'cell_b', 'cell_c', 'alpha', 'beta', 'gamma',
                 '#indexed_reflections', 'd_min_indexed',
                 'rmsd_x', 'rmsd_y', 'rmsd_phi')]
  for i in range(len(cell_params)):
    table_data.append((sweep_dir_cryst[i],
                       str(cell_params[i][0]),
                       str(cell_params[i][1]),
                       str(cell_params[i][2]),
                       str(cell_params[i][3]),
                       str(cell_params[i][4]),
                       str(cell_params[i][5]),
                       str(n_indexed[i]),
                       str(d_min_indexed[i]),
                       str(rmsds[i][0]),
                       str(rmsds[i][1]),
                       str(rmsds[i][2]),
                       ))

  with open('results_indexed.txt', 'wb') as f:
    print >> f, table_utils.format(
      table_data, has_header=True, justify='right')

  cell_a = flex.double([params[0] for params in cell_params])
  cell_b = flex.double([params[1] for params in cell_params])
  cell_c = flex.double([params[2] for params in cell_params])
  cell_alpha = flex.double([params[3] for params in cell_params])
  cell_beta = flex.double([params[4] for params in cell_params])
  cell_gamma = flex.double([params[5] for params in cell_params])

  from matplotlib import pyplot
  from matplotlib.backends.backend_pdf import PdfPages

  pyplot.rc('font', family='serif')
  pyplot.rc('font', serif='Times New Roman')

  red, blue = '#B2182B', '#2166AC'
  hist = flex.histogram(n_strong_spots_dmin_4.as_double(), n_slots=20)
  hist.show()
  fig = pyplot.figure()
  ax = fig.add_subplot(1,1,1)
  ax.bar(hist.slot_centers(), hist.slots(), width=0.75*hist.slot_width(),
         color=blue, edgecolor=blue)
  ax.set_xlabel('Spot count')
  ax.set_ylabel('Frequency')
  pdf = PdfPages("spot_count_histogram.pdf")
  pdf.savefig(fig)
  pdf.close()
  #pyplot.show()

  hist = flex.histogram(n_indexed_lattices.as_double(),
                        n_slots=flex.max(n_indexed_lattices))
  hist.show()
  fig = pyplot.figure()
  ax = fig.add_subplot(1,1,1)
  ax.bar(range(int(hist.data_max())), hist.slots(),
         width=0.75*hist.slot_width(), align='center',
         color=blue, edgecolor=blue)
  ax.set_xlim(-0.5, hist.data_max()-0.5)
  ax.set_xticks(range(0,int(hist.data_max())))
  ax.set_xlabel('Number of indexed lattices')
  ax.set_ylabel('Frequency')
  pdf = PdfPages("n_indexed_lattices_histogram.pdf")
  pdf.savefig(fig)
  pdf.close()
  #pyplot.show()

  if flex.max(n_integrated_lattices) > 0:
    hist = flex.histogram(n_integrated_lattices.as_double(),
                          n_slots=flex.max(n_integrated_lattices))
    hist.show()
    fig = pyplot.figure()
    ax = fig.add_subplot(1,1,1)
    ax.bar(range(int(hist.data_max())), hist.slots(),
           width=0.75*hist.slot_width(),
           align='center', color=blue, edgecolor=blue)
    ax.set_xlim(-0.5, hist.data_max()-0.5)
    ax.set_xticks(range(0,int(hist.data_max())))
    ax.set_xlabel('Number of integrated lattices')
    ax.set_ylabel('Frequency')
    pdf = PdfPages("n_integrated_lattices_histogram.pdf")
    pdf.savefig(fig)
    pdf.close()
    #pyplot.show()

  fig, axes = pyplot.subplots(nrows=2, ncols=3, squeeze=False)
  for i, cell_param in enumerate(
    (cell_a, cell_b, cell_c, cell_alpha, cell_beta, cell_gamma)):
    ax = axes.flat[i]
    flex.min_max_mean_double(cell_param).show()
    print flex.median(cell_param)
    hist = flex.histogram(cell_param, n_slots=20)
    hist.show()
    ax.bar(hist.slot_centers(), hist.slots(), width=0.75*hist.slot_width(),
           color=blue, edgecolor=blue)
    ax.set_xlabel('Cell parameter')
    ax.set_ylabel('Frequency')
  pyplot.tight_layout()
  pdf = PdfPages("cell_parameters.pdf")
  pdf.savefig(fig)
  pdf.close()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
