from __future__ import division
from dials.util.options import OptionParser
from libtbx.phil import parse
import dxtbx, os
from dials.algorithms.peak_finding.spotfinder_factory import SpotFinderFactory
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory
from dials.algorithms.refinement import RefinerFactory
from dials.algorithms.profile_model.factory import ProfileModelFactory
from dials.algorithms.integration.integrator import IntegratorFactory
from dials.array_family import flex


phil_scope = parse('''
  input {
    single_img = None
      .type = str
      .help = Path to input image
  }
  output_dir = .
    .type = str
    .help = "Directory where results will be deposited"

  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
  indexing {
    include scope dials.algorithms.indexing.indexer.master_phil_scope
  }
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
''', process_includes=True)

from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
def __stupid_but_swig_safe__deepcopy__(self, memo):
  pass
cbf_wrapper.__deepcopy__ = __stupid_but_swig_safe__deepcopy__

def process_reference(reference):
    ''' Load the reference spots. '''
    from dials.util.command_line import Command
    from dials.array_family import flex
    if reference is None:
      return None
    assert("miller_index" in reference)
    Command.start('Removing reference spots with invalid coordinates')
    mask = flex.bool([x == (0, 0, 0) for x in reference['xyzcal.mm']])
    reference.del_selected(mask)
    mask = flex.bool([h == (0, 0, 0) for h in reference['miller_index']])
    reference.del_selected(mask)
    Command.end('Removed reference spots with invalid coordinates, %d remaining' %
                len(reference))
    return reference

def run():
  parser = OptionParser(
    phil = phil_scope)

  params, options = parser.parse_args(show_diff_phil=True)
  assert params.input.single_img is not None
  assert params.output_dir is not None

  # load the image
  img = dxtbx.load(params.input.single_img)
  imgset = MemImageSet([img])
  datablock = DataBlockFactory.from_imageset(imgset)[0]

  spotfinder = SpotFinderFactory.from_parameters(params)
  reflections = spotfinder(datablock)

  base_name = os.path.splitext(params.input.single_img)[0]
  reflections.as_pickle(os.path.join(params.output_dir, base_name + "_strong.pickle"))

  # DGW commented out as reflections.minimum_number_of_reflections no longer exists
  #if len(reflections) < params.refinement.reflections.minimum_number_of_reflections:
  #  print "Not enough spots to index"
  #  return

  # create the spot finder

  print "Spotfinder spots found:", len(reflections)

  if params.indexing.method == "fft3d":
    from dials.algorithms.indexing.fft3d import indexer_fft3d as indexer
  elif params.indexing.method == "fft1d":
    from dials.algorithms.indexing.fft1d import indexer_fft1d as indexer
  elif params.method == "real_space_grid_search":
    from dials.algorithms.indexing.real_space_grid_search \
         import indexer_real_space_grid_search as indexer
  try:
    idxr = indexer(reflections, [imgset], params=params.indexing)
  except (RuntimeError, Sorry) as e:
    print str(e)
    return

  indexed = idxr.refined_reflections
  experiments = idxr.refined_experiments
  #from dxtbx.model.experiment.experiment_list import ExperimentListDumper
  #dump = ExperimentListDumper(experiments)
  #dump.as_json(os.path.join(params.output_dir, base_name + "_experiments.json"))
  indexed.as_pickle(os.path.join(params.output_dir, base_name + "_indexed.pickle"))

  refiner = RefinerFactory.from_parameters_data_experiments(
    params, indexed, experiments)

  refiner.run()
  refined_experiments = refiner.get_experiments()
  #dump = ExperimentListDumper(refined_experiments)
  #dump.as_json(os.path.join(params.output_dir, base_name + "_refined.json"))

  # Compute the profile model
  # Predict the reflections
  # Match the predictions with the reference
  # Create the integrator
  reference = indexed

  reference = process_reference(reference)
  profile_model = ProfileModelFactory.create(params, refined_experiments, reference)
  predicted = flex.reflection_table.from_predictions_multi(
    refined_experiments,
    dmin=params.prediction.dmin,
    dmax=params.prediction.dmax,
    margin=params.prediction.margin,
    force_static=params.prediction.force_static)
  predicted.match_with_reference(reference)
  integrator = IntegratorFactory.create(params, experiments, profile_model, predicted)

  # Integrate the reflections
  integrated = integrator.integrate()
  integrated.as_pickle(os.path.join(params.output_dir, base_name + "_integrated.pickle"))

if __name__ == "__main__":
  run()
