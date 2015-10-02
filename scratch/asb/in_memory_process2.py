from __future__ import division
from dials.util.options import OptionParser
from libtbx.phil import parse
import dxtbx, os
from dxtbx.imageset import MemImageSet
from dxtbx.datablock import DataBlockFactory

phil_scope = parse('''
  input {
    single_img = None
      .type = str
      .help = Path to input image
  }
  output {
    datablock_filename = %s_datablock.json
      .type = str
      .help = "The filename for output datablock"

    strong_filename = %s_strong.pickle
      .type = str
      .help = "The filename for strong reflections from spot finder output."

    indexed_filename = %s_indexed.pickle
      .type = str
      .help = "The filename for indexed reflections."

    refined_experiments_filename = %s_refined_experiments.json
      .type = str
      .help = "The filename for saving refined experimental models"

    integrated_filename = %s_integrated.pickle
      .type = str
      .help = "The filename for final integrated reflections."

    profile_filename = None
      .type = str
      .help = "The filename for output reflection profile parameters"
  }

  include scope dials.algorithms.peak_finding.spotfinder_factory.phil_scope
  indexing {
    include scope dials.algorithms.indexing.indexer.master_phil_scope
  }
  include scope dials.algorithms.refinement.refiner.phil_scope
  include scope dials.algorithms.integration.integrator.phil_scope
  include scope dials.algorithms.profile_model.factory.phil_scope
  include scope dials.algorithms.spot_prediction.reflection_predictor.phil_scope
''', process_includes=True)

from xfel.cftbx.detector.cspad_cbf_tbx import cbf_wrapper
def __stupid_but_swig_safe__deepcopy__(self, memo):
  pass
cbf_wrapper.__deepcopy__ = __stupid_but_swig_safe__deepcopy__

from xfel.command_line.dials_process import Script as DialsProcessScript
class InMemScript(DialsProcessScript):
  def __init__(self):
    self.parser = OptionParser(
      phil = phil_scope)

  def run(self):
    params, options = self.parser.parse_args(show_diff_phil=True)
    assert params.input.single_img is not None

    filebase = os.path.splitext(params.input.single_img)[0]

    for item in dir(params.output):
      value = getattr(params.output, item)
      try:
        if "%s" in value:
          setattr(params.output, item, value%filebase)
      except Exception:
        pass

    self.params = params
    self.options = options

    # load the image
    img = dxtbx.load(params.input.single_img)
    imgset = MemImageSet([img])
    datablock = DataBlockFactory.from_imageset(imgset)[0]

    # Cannot export MemImageSets
    #if self.params.output.datablock_filename:
      #from dxtbx.datablock import DataBlockDumper
      #dump = DataBlockDumper(datablock)
      #dump.as_json(self.params.output.datablock_filename)

    observed = self.find_spots(datablock)
    experiments, indexed = self.index(datablock, observed)
    experiments = self.refine(experiments, indexed)
    integrated = self.integrate(experiments, indexed)


if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = InMemScript()
    script.run()
  except Exception as e:
    halraiser(e)

