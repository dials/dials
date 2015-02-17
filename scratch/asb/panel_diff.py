from __future__ import division
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks
from dials.util.options import flatten_experiments
from scitbx.matrix import col

"""
Given two .json files with dxtbx detectors in them, show the difference
between them at each hierarchy level for the first panel.

Example usage: libtbx.python panel_diff.py json1 json2
"""

class Script(object):
  def __init__(self):
    # Create the parser
    self.parser = OptionParser(
      read_experiments=True,
      read_datablocks=True,
      read_reflections=True,
      read_datablocks_from_images=True,
      check_format=False)

  def run(self):
    # load at least two detectors from the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    datablocks  = flatten_datablocks(params.input.datablock)
    experiments = flatten_experiments(params.input.experiments)

    # collect all detectors found
    all_detectors = []
    for db in datablocks:
      all_detectors.extend(db.unique_detectors())

    all_detectors.extend(experiments.detectors())

    assert len(all_detectors) >= 2

    a = all_detectors[0]
    b = all_detectors[1]

    level = 0
    pga = a.hierarchy()
    pgb = b.hierarchy()

    while True:
      # starting at the top of the hierarchy, show diffs in local origins at each level
      print "Level", level
      oa = col(pga.get_local_origin())
      ob = col(pgb.get_local_origin())

      print "  Detector a", oa.elems
      print "  Detector b", ob.elems
      print "  Diffs", (ob-oa).elems

      if hasattr(pga, 'children'):
        pga = pga[0]
        pgb = pgb[0]
        level += 1
      else:
        break

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
