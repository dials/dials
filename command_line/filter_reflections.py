#!/usr/bin/env python
#
# dials.filter_reflections.py
#
#  Copyright (C) 2015 STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.filter_reflections

from __future__ import absolute_import, division
from dials.array_family import flex

help_message = '''

This program takes reflection files as input and filters them based on user-
specified criteria, to write out a subset of the original file.

Currently, only filtering by reflection flags is supported. Inclusions are
processed first and are combined by logical OR. Specifying no inclusions is a
special case in which all reflections are included for further filtering.
Exclusions are then processed, such that reflections are removed only if none of
the specified flags are set. Different results, such as filtering to include
only reflections with both flag1 AND flag2 set may be achieved by multiple runs
of the program.

Example::

  dev.dials.filter_reflections refined.pickle inclusions.flag=used_in_refinement

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

      inclusions {
        flag = %s
          .type = choice
          .help = "Include reflections with this flag to form the working set."
          .multiple = True
      }

      exclusions {
        flag = %s
          .type = choice
          .help = "Exclude reflections from the working set with this flag."
          .multiple = True
      }

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

    ''' % tuple([' '.join(self.flag_names)] * 2)

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
      print st.format()
    else:
      print "No flags set"

    return

  def run(self):
    '''Execute the script.'''
    from dials.array_family import flex
    from dials.util.options import flatten_reflections
    from dials.util.options import flatten_datablocks
    from dials.util.options import flatten_experiments
    from libtbx.utils import Sorry

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    reflections = flatten_reflections(params.input.reflections)

    if params.input.datablock is not None:
      datablocks = flatten_datablocks(params.input.datablock)
      assert len(datablocks) == 1
      imagesets = datablocks[0].extract_imagesets()
      assert len(imagesets) == 1
      imageset = imagesets[0]
    elif params.input.experiments is not None:
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

    print "{0} reflections loaded".format(len(reflections))

    if (len(params.inclusions.flag) == 0 and
        len(params.exclusions.flag) == 0 and
        params.d_min is None and params.d_max is None and
        params.partiality.min is None and params.partiality.max is None and
        not params.ice_rings.filter):
      print "No filter specified. Performing analysis instead."
      return self.analysis(reflections)

    # Build up the initial inclusion selection
    inc = flex.bool(len(reflections), True)
    # 2016/07/06 GW logic here not right should be && for each flag not or?
    for flag in params.inclusions.flag:
      sel = reflections.get_flags(getattr(reflections.flags, flag))
      inc = inc & sel
    reflections = reflections.select(inc)

    print "{0} reflections selected to form the working set".format(len(reflections))

    # Make requested exclusions from the current selection
    exc = flex.bool(len(reflections))
    for flag in params.exclusions.flag:
      print flag
      sel = reflections.get_flags(getattr(reflections.flags, flag))
      exc = exc | sel
    reflections = reflections.select(~exc)

    print "{0} reflections excluded from the working set".format(exc.count(True))

    # Filter based on resolution
    if params.d_min is not None:
      selection = reflections['d'] >= params.d_min
      reflections = reflections.select(selection)
      print "Selected %d reflections with d >= %f" % (len(reflections), params.d_min)

    # Filter based on resolution
    if params.d_max is not None:
      selection = reflections['d'] <= params.d_max
      reflections = reflections.select(selection)
      print "Selected %d reflections with d <= %f" % (len(reflections), params.d_max)

    # Filter based on partiality
    if params.partiality.min is not None:
      selection = reflections['partiality'] >= params.partiality.min
      reflections = reflections.select(selection)
      print "Selected %d reflections with partiality >= %f" % (len(reflections), params.partiality.min)

    # Filter based on partiality
    if params.partiality.max is not None:
      selection = reflections['partiality'] <= params.partiality.max
      reflections = reflections.select(selection)
      print "Selected %d reflections with partiality <= %f" % (len(reflections), params.partiality.max)

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

      print "Rejecting %i reflections at ice ring resolution" %ice_sel.count(True)
      reflections = reflections.select(~ice_sel)
      #reflections = reflections.select(ice_sel)

    # Save filtered reflections to file
    if params.output.reflections:
      print "Saving {0} reflections to {1}".format(len(reflections),
                                                   params.output.reflections)
      reflections.as_pickle(params.output.reflections)

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
