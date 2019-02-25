#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# LIBTBX_SET_DISPATCHER_NAME dev.dials.combine_datablocks

"""
Utility script to combine multiple datablock and strong spot files into
one datablock and one strong spot file. Imagesets are matched to reflections
in the order they are provided as input.

"""

#!/usr/bin/env dials.python
from __future__ import absolute_import, division, print_function

from dials.array_family import flex
from dials.util.options import flatten_datablocks, flatten_reflections
from dxtbx.datablock import DataBlock, DataBlockDumper
from libtbx.phil import parse
from dials.util import Sorry

# The phil scope
phil_scope = parse('''

  output {
    datablocks_filename = combined_datablocks.json
      .type = str
      .help = "The filename for combined datablocks"

    reflections_filename = combined_strong.pickle
      .type = str
      .help = "The filename for combined reflections"

    compact = False
      .type = bool
      .help = "For JSON output use compact representation"

  }
''', process_includes=True)

help_message = '''

Utility script to combine multiple strong spot and datablock files into
one multi-imageset strong spot file and one datablock file. Imagesets are
matched to reflections in the order they are provided as input.

Examples::

  dev.dials.combine_datablocks datablock_0.json datablock_1.json \\
    strong_0.pickle strong_1.pickle

'''

class Script(object):

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    # The script usage
    usage  = "usage: %s [options] [param.phil] " \
             "experiments1.json experiments2.json reflections1.pickle " \
             "reflections2.pickle..." \
             % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_datablocks=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)

    try:
      assert len(params.input.reflections) == len(params.input.datablock)
    except AssertionError:
      raise Sorry("The number of input reflections files does not match the "
        "number of input datablocks")

    datablocks = flatten_datablocks(params.input.datablock)
    reflections = flatten_reflections(params.input.reflections)

    if len(reflections):
      r = self.combine_reflections(reflections)
      # print number of reflections per imageset
      from libtbx.table_utils import simple_table
      max_id = max(r['id'])
      header = ["Imageset", "Nref"]
      nrefs_per_imset = [(r['id'] == i).count(True) for i in range(max_id + 1)]
      rows = [(str(i), str(n)) for (i, n) in enumerate(nrefs_per_imset)]
      st = simple_table(rows, header)
      print(st.format())
      rf = params.output.reflections_filename
      print('Saving combined reflections to {0}'.format(rf))
      r.as_pickle(rf)

    if len(datablocks):
      db = self.combine_datablocks(datablocks)
      dbf = params.output.datablocks_filename
      print('Saving combined datablocks to {0}'.format(dbf))
      dump = DataBlockDumper(db)
      dump.as_file(dbf, compact=params.output.compact)

    return

  def combine_datablocks(self, datablocks):
    imageset_list=[]
    for d in datablocks:
      imageset_list.extend(d.extract_imagesets())
    return DataBlock(imageset_list)

  def combine_reflections(self, reflections):
    new_reflections = flex.reflection_table()
    start_id = 0
    for rt in reflections:
      id_range = max(rt['id']) - min(rt['id']) + 1
      if len(set(rt['id'])) != id_range:
        raise Sorry('imageset ids not contiguous')
      rt['id'] = rt['id'] - min(rt['id']) + start_id
      start_id = max(rt['id']) + 1
      new_reflections.extend(rt)
    return new_reflections

if __name__ == "__main__":
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
