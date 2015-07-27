#!/usr/bin/env dials.python

"""This script is intended for quadrant refinement using still shot data
collected on a CSPAD detector. This version of the script uses a 'hybrid
minimiser'. Rather than a single joint refinement job of all crystals and the
detector, only the detector parameters are refined at first (using all data)
then each crystal is refined individually. This forms one macrocycle."""

from __future__ import division
import sys
sys.path.insert(1,'/reg/common/package/mpi4py/mpi4py-1.3.1/install/lib/python')
from mpi4py import MPI

from libtbx.phil import command_line, parse
from dxtbx.serialize import load as load_dxtbx
from dxtbx.model.experiment.experiment_list import ExperimentList, Experiment

from dials.model.serialize import load as load_dials
from dials.array_family import flex
from dials.algorithms.refinement import RefinerFactory

def load_input(exp_path, ref_path):

  refs = load_dials.reflections(ref_path)
  exp = load_dxtbx.experiment_list(exp_path , check_format=False)[0]
  return refs, exp

class ExperimentFromCrystal(object):

  def __init__(self, reference_beam, reference_detector):

    self.reference_beam = reference_beam
    self.reference_detector = reference_detector
    return

  def __call__(self, crystal):

    return Experiment(beam=self.reference_beam,
                      detector=self.reference_detector,
                      crystal=crystal)

class DetectorRefiner(object):

  user_phil=parse("""
  refinement{
    parameterisation {
      beam.fix=all
      crystal.fix=all
      detector.hierarchy_level=1
      sparse=True
    }
    target.gradient_calculation_blocksize=100000
    reflections.outlier.algorithm=tukey
    refinery.engine=LBFGScurvs
  }""")
  from dials.algorithms.refinement.refiner import phil_scope as refinement_phil
  working_phil = refinement_phil.fetch(sources=[user_phil])

  def __call__(self, experiments, reflections):

    self.working_phil.show()
    params = self.working_phil.extract()
    refiner = RefinerFactory.from_parameters_data_experiments(
        params, reflections, experiments, verbosity=1)
    refiner.run()
    experiments = refiner.get_experiments()
    print "Detector refinement finished"

    return experiments

class CrystalRefiners(object):
  user_phil=parse("""
  refinement{
    parameterisation {
      beam.fix=all
      detector.fix=all
    }
    reflections.outlier.algorithm=tukey
    refinery.engine=LBFGScurvs
  }""")
  from dials.algorithms.refinement.refiner import phil_scope as refinement_phil
  working_phil = refinement_phil.fetch(sources=[user_phil])

  def __call__(self, experiments, reflections):
    comm = MPI.COMM_WORLD
    rank = comm.Get_rank()
    params = self.working_phil.extract()

    if rank == 0:
      data = []
      size =comm.Get_size()
      chunk_size = len(experiments) // size
      remainder = len(experiments) % size
      pointer = 0

      self.working_phil.show()

      for i in xrange(size):
        if i < remainder:
          sel_range = xrange(pointer,pointer+chunk_size+1)
        else:
          sel_range = xrange(pointer,pointer+chunk_size)

        sel = flex.bool(len(reflections))
        for exp_id in sel_range:
          sel |= reflections['id'] == exp_id

        if i < remainder:
          data.append((range(pointer,pointer+chunk_size+1),experiments[pointer:pointer+chunk_size+1],reflections.select(sel)))
          pointer += 1
        else:
          data.append((range(pointer,pointer+chunk_size),experiments[pointer:pointer+chunk_size],reflections.select(sel)))
        pointer += chunk_size

    else:
      data = None

    data = comm.scatter(data, root=0)

    for i, (iexp, exp) in enumerate(zip(data[0],data[1])):

      print "Refining crystal", iexp
      # reflection subset for a single experiment
      refs = data[2].select(data[2]['id'] == iexp)
      refs['id'] = flex.size_t(len(refs),0)
      # experiment list for a single experiment
      exps=ExperimentList()
      exps.append(exp)
      refiner = RefinerFactory.from_parameters_data_experiments(
        params, refs, exps, verbosity=1)
      # do refinement
      refiner.run()
      refined_exps = refiner.get_experiments()
      # replace this experiment with the refined one
      data[1][i] = refined_exps[0]

    data = comm.gather(data, root=0)
    if rank == 0:
      for chunk in data:
        for iexp, experiment in zip(chunk[0], chunk[1]):
          experiments[iexp] = experiment

      return experiments
    else:
      assert data == None

if __name__ =="__main__":

  comm = MPI.COMM_WORLD
  rank = comm.Get_rank()


  if len(sys.argv) != 2: exit("please pass the path to a phil file")
  phil = sys.argv[1]

  master_phil = parse("""
    input
      .multiple = True
    {
      experiments = None
        .type = path
      reflections = None
        .type = path
    }
    n_macrocycles = 1
      .type = int(value_min=1)
    """)

  cmd_line = command_line.argument_interpreter(master_params=master_phil)
  working_phil = cmd_line.process_and_fetch(args=(phil,))
  working_params = working_phil.extract()

  assert len(working_params.input) > 1

  if rank == 0:
    for input in working_params.input:
      print input.experiments, input.reflections

    print len(working_params.input), "datasets specified as input"

    e = enumerate(working_params.input)
    i, line = e.next()
    reflections, exp = load_input(line.experiments, line.reflections)
    assert reflections['id'].all_eq(0)
    from dials.algorithms.indexing.indexer import indexer_base
    reflections = indexer_base.map_spots_pixel_to_mm_rad(reflections, exp.detector, exp.scan)
    experiment_from_crystal=ExperimentFromCrystal(exp.beam, exp.detector)

    experiments=ExperimentList()
    experiments.append(experiment_from_crystal(exp.crystal))

    for i, line in e:
      refs, exp = load_input(line.experiments, line.reflections)
      print i, line.reflections, len(refs)
      refs['id'] = flex.size_t(len(refs),i)
      refs = indexer_base.map_spots_pixel_to_mm_rad(refs, exp.detector, exp.scan)
      reflections.extend(refs)
      experiments.append(experiment_from_crystal(exp.crystal))

    dr = DetectorRefiner()
  else:
    experiments = None
    reflections = None
    dr = None
  cr = CrystalRefiners()


  for cycle in range(working_params.n_macrocycles):

    if rank == 0:
      print "MACROCYCLE %02d" % (cycle + 1)
      print "=============\n"
      # first run: multi experiment joint refinement of detector with fixed beam and
      # crystals
      experiments = dr(experiments, reflections)
    else:
      experiments = None

    # second run
    experiments = cr(experiments, reflections)

  if rank == 0:
    # save the refined experiments
    from dxtbx.model.experiment.experiment_list import ExperimentListDumper
    dump = ExperimentListDumper(experiments)
    experiments_filename = "refined_experiments.json"
    dump.as_json(experiments_filename)
    print "refined geometry written to {0}".format(experiments_filename)
