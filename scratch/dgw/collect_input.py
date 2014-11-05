#!/usr/bin/env python
#
#  Copyright (C) 2013 Diamond Light Source and STFC Rutherford Appleton
#  Laboratory, UK.
#
#  Author: David Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

"""Prepare input.phil for multi-experiment refinement by walking the filesystem
hierarchy rooted at the current working directory looking for files that
match the supplied template"""

from __future__ import division

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from dials.util.options import OptionParser
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''

      input {
        experiments_template = None
          .type = path
          .multiple = True

        reflections_template = None
          .type = path
          .multiple = True
      }

      output {
        phil_filename = experiments_and_reflections.phil
          .type = str
          .help = "The filename for combined experiments and reflections"
      }
    ''', process_includes=True)

    self.output_phil = parse('''
        input {
          experiments = None
            .type = path
            .multiple = True
            .help = "The experiment list file path"
          reflections = None
            .type = path
            .multiple = True
            .help = "The reflection table file path"
        }
      ''')

    # The script usage
    usage  = "usage: %s [options] [param.phil] " % libtbx.env.dispatcher_name

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope)

    return

  @staticmethod
  def find_file_pairs(e_templ, r_templ):
    import os

    pairs = []
    for dirpath, _, filenames in os.walk("."):
      dirpath = os.path.abspath(dirpath)
      for e, r in zip(e_templ, r_templ):
        if e in filenames and r in filenames:
          e_path = os.path.join(dirpath, e)
          r_path = os.path.join(dirpath, r)
          pairs.append((e_path, r_path))

    return pairs


  def run(self):
    '''Execute the script.'''

    # Parse the command line
    params, options = self.parser.parse_args(show_diff_phil=True)
    e_templ = params.input.experiments_template
    r_templ = params.input.reflections_template

    # Extract output_phil to a python object
    working = self.output_phil.extract()

    # remove first elts, which are None
    working.input.experiments = []
    working.input.reflections = []

    # do stuff
    pairs = self.find_file_pairs(e_templ, r_templ)

    for (e, r) in pairs:
      working.input.experiments.append(e)
      working.input.reflections.append(r)

    # display
    output_phil = self.output_phil.format(python_object=working)
    output_phil.show()

    # write
    with open(params.output.phil_filename, "w") as f:
      f.write(output_phil.as_str())

    return

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
