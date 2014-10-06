from __future__ import division

# Script to convert the output from refine_quadrants to a header file
# Note hardcoded distance of 100 isn't relevant for just a cbf header

from dials.util.options import OptionParser
from dials.util.options import flatten_experiments
from xfel.cftbx.detector.cspad_cbf_tbx import write_cspad_cbf, map_detector_to_basis_dict

class Script(object):
  def __init__(self):
    # Create the parser
    self.parser = OptionParser(
      read_experiments=True)

  def run(self):
    params, options = self.parser.parse_args(show_diff_phil=True)
    experiments = flatten_experiments(params.input.experiments)

    detector = experiments[0].detector

    metro = map_detector_to_basis_dict(detector)
    write_cspad_cbf(None, metro, 'cbf', None, 'refined_detector.def', None, 100, header_only=True)

    print "Done"

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

