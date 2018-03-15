from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

import iotbx.phil

from scitbx.array_family import flex
from dials.algorithms.symmetry.cosym.generate_test_data import generate_test_data
from dials.algorithms.symmetry.cosym import analyse_datasets


phil_scope = iotbx.phil.parse('''\

d_min = 1
  .type = float(value_min=0)

sigma = 0.1
  .type = float(value_min=0.0)
  .help = "Width of normal distribution used for adding random noise to intensities."

sample_size = 100
  .type = int(value_min=1)
  .help = "Number of datasets to generate."

unit_cell = None
  .type = unit_cell

unit_cell_volume = 1000
  .type = float(value_min=100)
  .help = "Volume of unit cell used to generate a compatible crystal symmetry"
          "if a specific unit cell is not provided."

space_group = P2
  .type = space_group

map_to_p1 = False
  .type = bool

include scope dials.algorithms.symmetry.cosym.phil_scope

seed = 0
  .type = int(value_min=0)

twin_fractions = None
  .type = floats(value_min=0, value_max=0.5)

''', process_includes=True)


def run(args):

  import random

  interp = phil_scope.command_line_argument_interpreter()
  params, unhandled = interp.process_and_fetch(
    args, custom_processor='collect_remaining')
  params = params.extract()

  if params.seed is not None:
    flex.set_random_seed(params.seed)
    random.seed(params.seed)

  from dials.util import log
  # Configure the logging
  log.config()

  datasets, expected_reindexing_ops = generate_test_data(
    space_group=params.space_group.group(),
    lattice_group=params.lattice_group,
    unit_cell=params.unit_cell,
    unit_cell_volume=params.unit_cell_volume,
    seed=params.seed,
    d_min=params.d_min,
    sigma=params.sigma,
    sample_size=params.sample_size,
    map_to_p1=params.map_to_p1,
    twin_fractions=params.twin_fractions)

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


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
