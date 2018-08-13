#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dev.dials.filter_reflections

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

  dev.dials.filter_reflections integrated.pickle \
    flag_expression=used_in_refinement

  dev.dials.filter_reflections integrated.pickle \
    flag_expression="integrated & ~reference_spot"

  dev.dials.filter_reflections integrated.pickle \
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

    if params.input.datablock is not None and len(params.input.datablock):
      datablocks = flatten_datablocks(params.input.datablock)
      assert len(datablocks) == 1
      imagesets = datablocks[0].extract_imagesets()
      assert len(imagesets) == 1
      imageset = imagesets[0]
    elif params.input.experiments is not None and len(params.input.experiments):
      experiments = flatten_experiments(params.input.experiments)
      assert len(datablocks) == 1
      imageset = experiments[0].imageset
    else:
      imageset = None

    if len(reflections) == 0:
      self.parser.print_help()
      raise Sorry('No valid reflection file given')
    if len(reflections) != 1:
      self.parser.print_help()
      raise Sorry('Exactly 1 reflection file must be specified')
    reflections = reflections[0]

    # Check params
    if params.d_min is not None and params.d_max is not None:
      if params.d_min > params.d_max:
        raise Sorry("d_min must be less than d_max")
    if params.d_min is not None or params.d_max is not None:
      if 'd' not in reflections:
        raise Sorry("Reflection table has no resolution information")

    # Check params
    if params.partiality.min is not None and params.partiality.max is not None:
      if params.min > params.max:
        raise Sorry("partiality.min must be less than partiality.d_max")
    if params.partiality.min is not None or params.partiality.max is not None:
      if 'partiality' not in reflections:
        raise Sorry("Reflection table has no partiality information")

    print("{0} reflections loaded".format(len(reflections)))

    # Check if any filter has been set using diff_phil
    filter_def = [o for o in self.parser.diff_phil.objects
                  if o.name not in ['input', 'output']]
    if not filter_def:
      print("No filter specified. Performing analysis instead.")
      return self.analysis(reflections)

    # Filter by logical expression using flags
    if params.flag_expression is not None:
      inc = self.eval_flag_expression(params.flag_expression, reflections)
      reflections = reflections.select(inc)

    print("Selected {0} reflections by flags".format(len(reflections)))

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

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
