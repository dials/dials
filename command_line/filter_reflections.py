#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dials.filter_reflections

from __future__ import absolute_import, division, print_function

from libtbx.utils import Sorry
from dials.array_family import flex

help_message = '''

This program takes reflection files as input and filters them based on user-
specified criteria, to write out a subset of the original file.

Filtering is first done by evaluating the optional boolean 'flag_expression'
using reflection flag values. The operators allowed are '&' for 'and', '|' for
'or', and '~' for 'not'. Expressions may contain nested sub-expressions using
parentheses.

Following this, optional additional filters are applied according to values in
the reflection table, such as by resolution or user-defined masks.

If a reflection file is passed in to the program but no filtering parameters
are set, a table will be printed, giving the flag values present in the
reflection file.


Examples::

  dials.filter_reflections integrated.pickle \
    flag_expression=used_in_refinement

  dials.filter_reflections integrated.pickle \
    flag_expression="integrated & ~reference_spot"

  dials.filter_reflections integrated.pickle \
    flag_expression="indexed & (failed_during_summation | failed_during_profile_fitting)"

'''

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env
    from operator import itemgetter

    flags = list(flex.reflection_table.flags.names.iteritems())
    flags.sort(key=itemgetter(0))
    self.flag_names, self.flag_values = zip(*flags)
    phil_str = '''

      output {
        reflections = 'filtered.pickle'
          .type = str
          .help = "The filtered reflections output filename"
      }

      flag_expression = None
        .type = str
        .help = "Boolean expression to select reflections based on flag values"

      id = None
        .type = ints(value_min=0)
        .help = "Select reflections by experiment IDs"

      panel = None
        .type = ints(value_min=0)
        .help = "Select reflections by panels they intersect"

      d_min = None
        .type = float
        .help = "The maximum resolution"

      d_max = None
        .type = float
        .help = "The minimum resolution"

      partiality {
        min = None
          .type = float(value_min = 0, value_max = 1)
          .help = "The minimum reflection partiality for inclusion."
        max = None
          .type = float(value_min = 0, value_max = 1)
          .help = "The maximum reflection partiality for inclusion."
      }

      select_good_intensities = False
        .type = bool
        .help = "Combined filter to select only fully integrated and"
                "trustworthy intensities"

      dead_time
      {
        value = 0
          .help = "Detector dead time in ms, assumed to be at the end of the"
                  "exposure time."
          .type = float(value_min=0)

        reject_fraction = 0
          .help = "Reject reflections which overlap by more than the given"
                  "fraction with the dead region of the image."
          .type = float(value_min=0, value_max=1)
      }

      include scope dials.util.masking.ice_rings_phil_scope

    '''

    phil_scope = parse(phil_str, process_includes=True)

    # The script usage
    usage  = "usage: %s [options] experiment.json" % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=help_message,
      read_reflections=True,
      read_experiments=True,
      read_datablocks=True)

  def analysis(self, reflections):
    '''Print a table of flags present in the reflections file'''

    from libtbx.table_utils import simple_table
    header = ['flag','nref']
    rows = []
    for name, val in zip(self.flag_names, self.flag_values):
      n = (reflections.get_flags(val)).count(True)
      if n > 0: rows.append([name, "%d" % n])
    if len(rows) > 0:
      st = simple_table(rows, header)
      print(st.format())
    else:
      print("No flags set")

    return

  def eval_flag_expression(self, expression, reflections):
    """Test a Boolean expression of reflection flags for validity then
    evaluate it"""

    import token
    from tokenize import generate_tokens, TokenError, untokenize
    from StringIO import StringIO

    result = []
    g = generate_tokens(StringIO(expression).readline)

    flags = list(flex.reflection_table.flags.names.iteritems())

    # define shorthand function
    def get_flag(flag):
      return reflections.get_flags(getattr(reflections.flags, flag))

    while True:

      # Extract next token, catching unmatched brackets
      try:
        toknum, tokval, _, _, _ = g.next()
      except TokenError:
        raise Sorry("errors found in {0}".format(expression))
      except StopIteration:
        break

      # Catch unwanted token types
      if toknum not in [token.OP, token.NAME, token.ENDMARKER]:
        raise Sorry("invalid tokens found in {0}".format(expression))

      # Catch unwanted operators
      if toknum is token.OP and tokval not in "()|&~":
        raise Sorry("unrecognised operators found in {0}".format(expression))

      # Catch unrecognised flag names
      if toknum is token.NAME and tokval not in self.flag_names:
        raise Sorry("unrecognised flag name: {0}".format(tokval))

      # Replace names with valid lookups in the reflection table
      if toknum is token.NAME:
        ('NAME', 'get_flags')
        ('OP', '(')
        ('STRING', "'indexed'")
        ('OP', ')')
        result.extend([(token.NAME, 'get_flag'),
                       (token.OP, '('),
                       (token.STRING, repr(tokval)),
                       (token.OP, ')')])
      else:
        result.append((toknum, tokval))

    # Evaluate and return the result
    return eval(untokenize(result), {'get_flag': get_flag, '__builtins__':None},{})

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_datablocks
    from dials.util.options import flatten_experiments

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)

    if len(reflections) == 0:
      self.parser.print_help()
      raise Sorry('No valid reflection file given')
    if len(reflections) != 1:
      self.parser.print_help()
      raise Sorry('Exactly 1 reflection file must be specified')
    reflections = reflections[0]

    experiments = flatten_experiments(params.input.experiments)
    datablocks = flatten_datablocks(params.input.datablock)
    if len(experiments) > 1 or datablocks:
      if [len(datablocks), len(experiments)].count(1) != 1:
          self.parser.print_help()
          raise Sorry("Either a datablock or an experiment list may be provided"
                      " but not both together.")
    if datablocks:
      datablock = datablocks[0]
      imagesets = datablock.extract_imagesets()
    if len(experiments) > 1:
      imagesets = experiments.imagesets()

    # Check if any filter has been set using diff_phil
    filter_def = [o for o in self.parser.diff_phil.objects
                  if o.name not in ['input', 'output']]
    if not filter_def:
      print("No filter specified. Performing analysis instead.")
      return self.analysis(reflections)

    # Check params
    if params.d_min is not None and params.d_max is not None:
      if params.d_min > params.d_max:
        raise Sorry("d_min must be less than d_max")
    if params.d_min is not None or params.d_max is not None:
      if 'd' not in reflections:
        if len(experiments) > 0:
          print("Reflection table does not have resolution information. "
                "Attempting to calculate this from the experiment list")
          sel = reflections['id'] >= 0
          if sel.count(False) > 0:
            print("Removing {0} reflections with negative experiment id".format(
                sel.count(False)))
          reflections = reflections.select(sel)
          reflections.compute_d(experiments)
        else:
          raise Sorry("reflection table has no resolution information "
                      "and no experiment list provided to calculate it")

    # Check params
    if params.partiality.min is not None and params.partiality.max is not None:
      if params.min > params.max:
        raise Sorry("partiality.min must be less than partiality.d_max")
    if params.partiality.min is not None or params.partiality.max is not None:
      if 'partiality' not in reflections:
        raise Sorry("Reflection table has no partiality information")

    print("{0} reflections loaded".format(len(reflections)))

    # Filter by logical expression using flags
    if params.flag_expression is not None:
      inc = self.eval_flag_expression(params.flag_expression, reflections)
      reflections = reflections.select(inc)

    print("Selected {0} reflections by flags".format(len(reflections)))

    # Filter based on experiment ID
    if params.id:
      selection = reflections['id'] == params.id[0]
      for exp_id in params.id[1:]:
        selection = selection | (reflections['id'] == exp_id)
      reflections = reflections.select(selection)
      print("Selected %d reflections by experiment id" % (len(reflections)))

    # Filter based on panel number
    if params.panel:
      selection = reflections['panel'] == params.panel[0]
      for pnl_id in params.panel[1:]:
        selection = selection | (reflections['panel'] == pnl_id)
      reflections = reflections.select(selection)
      print("Selected %d reflections by panel number" % (len(reflections)))

    # Filter based on resolution
    if params.d_min is not None:
      selection = reflections['d'] >= params.d_min
      reflections = reflections.select(selection)
      print("Selected %d reflections with d >= %f" % (len(reflections), params.d_min))

    # Filter based on resolution
    if params.d_max is not None:
      selection = reflections['d'] <= params.d_max
      reflections = reflections.select(selection)
      print("Selected %d reflections with d <= %f" % (len(reflections), params.d_max))

    # Filter based on partiality
    if params.partiality.min is not None:
      selection = reflections['partiality'] >= params.partiality.min
      reflections = reflections.select(selection)
      print("Selected %d reflections with partiality >= %f" % (len(reflections), params.partiality.min))

    # Filter based on partiality
    if params.partiality.max is not None:
      selection = reflections['partiality'] <= params.partiality.max
      reflections = reflections.select(selection)
      print("Selected %d reflections with partiality <= %f" % (len(reflections), params.partiality.max))

    # 'Good' intensity selection
    if params.select_good_intensities:
      reflections = select_good_intensities(reflections)

    # Dead time filter
    if params.dead_time.value > 0:
      reflections = filter_by_dead_time(reflections, experiments,
          params.dead_time.value, params.dead_time.reject_fraction)

    # Filter powder rings
    if params.ice_rings.filter:
      from dials.algorithms.integration import filtering
      if 'd' in reflections:
        d_spacings = reflections['d']
      else:
        from cctbx import uctbx
        if 'rlp' not in reflections:
          assert imageset is not None
          from dials.algorithms.spot_finding.per_image_analysis import map_to_reciprocal_space
          reflections = map_to_reciprocal_space(reflections, imageset)
        d_star_sq = flex.pow2(reflections['rlp'].norms())
        d_spacings = uctbx.d_star_sq_as_d(d_star_sq)

      d_min = params.ice_rings.d_min
      width = params.ice_rings.width

      if d_min is None:
        d_min = flex.min(d_spacings)

      ice_filter = filtering.PowderRingFilter(
        params.ice_rings.unit_cell, params.ice_rings.space_group.group(), d_min, width)

      ice_sel = ice_filter(d_spacings)

      print("Rejecting %i reflections at ice ring resolution" %ice_sel.count(True))
      reflections = reflections.select(~ice_sel)

    # Save filtered reflections to file
    if params.output.reflections:
      print("Saving {0} reflections to {1}".format(len(reflections),
                                                   params.output.reflections))
      reflections.as_pickle(params.output.reflections)

    return

def select_good_intensities(integrated_data):
  """Combined filter, originally in dev.dials.filter_good_intensities"""
  if not (integrated_data.has_key('id') and
          integrated_data.has_key('intensity.sum.variance')):
    raise Sorry("reflection file is missing required keys")

  if not (min(integrated_data['id']) == max(integrated_data['id']) == 0):
    raise Sorry("only a single experiment is supported in this mode")

  selection = integrated_data['intensity.sum.variance'] <= 0
  if selection.count(True) > 0:
    integrated_data.del_selected(selection)
    print('Removing %d reflections with negative variance' % \
          selection.count(True))

  if 'intensity.prf.variance' in integrated_data:
    selection = integrated_data['intensity.prf.variance'] <= 0
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      print('Removing %d profile reflections with negative variance' % \
            selection.count(True))

  if 'partiality' in integrated_data:
    selection = integrated_data['partiality'] < 0.99
    if selection.count(True) > 0:
      integrated_data.del_selected(selection)
      print('Removing %d incomplete reflections' % \
        selection.count(True))

  return integrated_data

def filter_by_dead_time(reflections, experiments,
    dead_time=0, reject_fraction=0):
  """Combined filter, originally in dev.dials.filter_dead_time"""

  if len(experiments) == 0:
    raise Sorry("an experiment list must be provided to filter by dead time")

  if len(experiments) > 1:
    raise Sorry("only a single experiment is supported in this mode")
  experiment = experiments[0]

  sel = reflections.get_flags(reflections.flags.integrated)
  reflections = reflections.select(sel)

  if len(reflections) == 0:
    raise Sorry("no integrated reflections present")

  goniometer = experiment.goniometer
  beam = experiment.beam

  m2 = goniometer.get_rotation_axis()
  s0 = beam.get_s0()

  from dials.array_family import flex
  phi1 = flex.double()
  phi2 = flex.double()

  phi_range = reflections.compute_phi_range(
    goniometer.get_rotation_axis(),
    beam.get_s0(),
    experiment.profile.sigma_m(deg=False),
    experiment.profile.n_sigma())
  phi1, phi2 = phi_range.parts()

  scan = experiment.scan
  exposure_time = scan.get_exposure_times()[0]
  assert scan.get_exposure_times().all_eq(exposure_time)
  phi_start, phi_width = scan.get_oscillation(deg=False)
  phi_range_dead = phi_width * (dead_time/1000) / exposure_time

  sel_good = flex.bool(len(reflections), True)

  start, end = scan.get_array_range()
  for i in range(start, end):
    phi_dead_start = phi_start + (i+1) * phi_width - phi_range_dead
    phi_dead_end = phi_dead_start + phi_range_dead

    left = phi1.deep_copy()
    left.set_selected(left < phi_dead_start, phi_dead_start)

    right = phi2.deep_copy()
    right.set_selected(right > phi_dead_end, phi_dead_end)

    overlap = (right - left)/(phi2-phi1)

    sel = overlap > reject_fraction

    sel_good.set_selected(sel, False)
    print('Rejecting %i reflections from image %i' %(sel.count(True), i))

  print('Keeping %i reflections (rejected %i)' %(
    sel_good.count(True), sel_good.count(False)))

  return reflections.select(sel_good)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
