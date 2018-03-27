#!/usr/bin/env dials.python
from __future__ import absolute_import, division
from __future__ import print_function
from libtbx.phil import parse
from scitbx import matrix
from libtbx.table_utils import simple_table

# LIBTBX_SET_DISPATCHER_NAME dev.dials.compare_beam_spindle_angles

help_message = """

A script to calculate the angle between the beam and the rotation axis for
multiple experiments, and report the mean and standard deviation of these.
This is a utility for investigating mechanical properties of various goniometer
systems used at Diamond Light Source and is a partial response to
https://github.com/dials/dials/issues/271

Example::

  dials.compare_beam_spindle_angles experiments1.json experiments2.json ...

"""

class Script(object):

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''

    print_precision = 4
      .type=int(value_min=0)
      .help="Number of decimal places to print values with"

    ''', process_includes=True)

    # The script usage
    usage  = ("usage: %s [options] [param.phil] experiments1.json "
              "experiments2.json..." % libtbx.env.dispatcher_name)

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    '''Execute the script.'''

    from dials.util.options import flatten_experiments
    from libtbx.utils import Sorry
    from dials.array_family import flex

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    # Try to load the experiments
    if not params.input.experiments:
      print("No Experiments found in the input")
      self.parser.print_help()
      return

    experiments = flatten_experiments(params.input.experiments)
    print("{0} experiments loaded".format(len(experiments)))

    us0_vecs = self.extract_us0_vecs(experiments)
    e_vecs = self.extract_rotation_axes(experiments)

    angles = [us0.angle(e, deg=True) for us0, e in zip(us0_vecs, e_vecs)]

    fmt = "{:." + str(params.print_precision) + "f}"
    header = ['Exp\nid','Beam direction', 'Rotation axis', 'Angle (deg)']
    rows = []
    for iexp, (us0, e, ang) in enumerate(zip(us0_vecs, e_vecs, angles)):
      beam_str = " ".join([fmt] * 3).format(*us0.elems)
      e_str = " ".join([fmt] * 3).format(*e.elems)
      rows.append([str(iexp), beam_str, e_str, fmt.format(ang)])
    if len(rows) > 0:
      st = simple_table(rows, header)
      print(st.format())

    # mean and sd
    if len(rows) > 1:
      angles = flex.double(angles)
      mv = flex.mean_and_variance(angles)

      print("Mean and standard deviation of the angle")
      print (fmt.format(mv.mean()) + " +/- " + fmt.format(
        mv.unweighted_sample_standard_deviation()))
      print()

    return

  def extract_us0_vecs(self, experiments):
    return [matrix.col(e.beam.get_unit_s0()) for e in experiments]

  def extract_rotation_axes(self, experiments):
    axes = []
    for iexp, exp in enumerate(experiments):
      try:
        axes.append(matrix.col(exp.goniometer.get_rotation_axis()))
      except AttributeError:
        raise Sorry("Experiment with id {0} has no goniometer".format(iexp))
    return axes

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)

