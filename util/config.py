#!/usr/bin/env python
#
# config.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

class CompletionGenerator(object):
  '''A class to auto-generate bash-completion hints for dials'''

  def __init__(self):
    '''Initialise the class with the list of programs'''
    import libtbx.load_env
    import os

    # Find the dials source directory
    self.dist_path = libtbx.env.dist_path('dials')

    # Find the dials distribution directory
    build_path = abs(libtbx.env.build_path)

    # Set the location of the output directory
    self.output_directory = os.path.join(
        build_path, 'dials', 'autocomplete')
    try:
      os.makedirs(self.output_directory)
    except OSError:
      pass

  def generate(self):
    '''Generate the autocompletion hints.'''
    import os
    commands_dir = os.path.join(self.dist_path, 'command_line')
    print 'Generating command line completion hints for',
    for file in sorted(os.listdir(commands_dir)):
      if not file.startswith('_') and file.endswith('.py'):
        if 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in open(os.path.join(commands_dir, file)).read():
          command_name = 'dials.%s' % file[:-3]
          print command_name,
          self._generate_single(commands_dir, file, command_name)
    print

  def _generate_single(self, directory, executable, command):
    '''Generate a hints file for a single program.'''
    import os
    import subprocess

    # Save the generated hints to file
    with open(os.path.join(self.output_directory, command), 'w') as output:
      returncode = subprocess.call(["libtbx.python", os.path.join(directory, executable), "--export-autocomplete-hints"], stdout=output)
      print ("[OK]" if returncode == 0 else "[FAIL:%d]" % returncode),

if __name__ == '__main__':
  gen = CompletionGenerator()
  gen.generate()
