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
      for file in os.listdir(self.output_directory):
        file_path = os.path.join(self.output_directory, file)
        try:
          if os.path.isfile(file_path):
            os.unlink(file_path)
        except Exception, e:
          pass
      pass

  def generate(self):
    '''Generate the autocompletion init script.'''
    command_list = []

    import os
    commands_dir = os.path.join(self.dist_path, 'command_line')
    print 'Identifying autocompletable commands:',
    for file in sorted(os.listdir(commands_dir)):
      if not file.startswith('_') and file.endswith('.py'):
        if 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in open(os.path.join(commands_dir, file)).read():
          command_name = 'dials.%s' % file[:-3]
          print command_name,
          command_list.append(command_name)
    print
    with open(os.path.join(self.output_directory, 'init.sh'), 'w') as loader:
      loader.write('''#!/bin/bash
if [ ! -e "%s" ]; then touch "%s" 2>/dev/null; if [ -e "%s" ]; then
for cmd in %s; do
 echo Generating command line completion hints for $cmd
 $cmd --export-autocomplete-hints > "%s" || rm "%s"
done; fi; fi
source %s/util/autocomplete.sh
''' % (
        os.path.join(self.output_directory, 'runonce'),
        os.path.join(self.output_directory, 'runonce'),
        os.path.join(self.output_directory, 'runonce'),
        " ".join(command_list),
        os.path.join(self.output_directory, '${cmd}'),
        os.path.join(self.output_directory, '${cmd}'),
        self.dist_path))

if __name__ == '__main__':
  gen = CompletionGenerator()
  gen.generate()
