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

    # Set the location of the output directory
    self.output_directory = libtbx.env.under_build(os.path.join('dials', 'autocomplete'))
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
    import os
    commands_dir = os.path.join(self.dist_path, 'command_line')
    command_list = []
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
for cmd in %s; do
 if [ ! -e "%s" ]; then
  echo Generating command line completion hints for $cmd
  $cmd --export-autocomplete-hints > "%s" || rm "%s"
 fi
done
source %s/util/autocomplete.sh
%s
''' % (
        " ".join(command_list),
        os.path.join(self.output_directory, '${cmd}'),
        os.path.join(self.output_directory, '${cmd}'),
        os.path.join(self.output_directory, '${cmd}'),
        self.dist_path,
        "\n".join(["complete -F _dials_autocomplete %s" % cmd for cmd in command_list])
     ))
    with open(os.path.join(self.output_directory, 'SConscript'), 'w') as builder:
      builder.write('''
import os.path
import libtbx.load_env
def dispatcher_outer(name):
  return os.path.join(libtbx.env.under_build('bin'), name)
def dispatcher_inner(name):
  return os.path.join(libtbx.env.dist_path('dials'), 'command_line', '%%s.py' %% name.partition('.')[2])
env = Environment()
env.Append( BUILDERS={'AutoComplete': Builder(action='$SOURCE --export-autocomplete-hints > $TARGET')} )
for cmd in [%s]:
  env.AutoComplete(cmd, [dispatcher_outer(cmd), dispatcher_inner(cmd)])
''' % ( ', '.join(["'%s'" % cmd for cmd in command_list]) ))

  def install(self):
    '''Permanently install the autocompletion init script into setpaths-scripts.'''
    import libtbx.load_env
    import os
    build_path = abs(libtbx.env.build_path)
    print "Installing autocompletion script into:",
    for file in os.listdir(build_path):
      if file.startswith('setpath') and file.endswith('.sh'):
        if not 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in open(os.path.join(build_path, file)).read():
          print file,
          with open(os.path.join(build_path, file), 'a') as script:
            script.write('\n\n# DIALS_ENABLE_COMMAND_LINE_COMPLETION\n')
            script.write('[ -z "$BASH_VERSIONINFO" ] && source %s\n' % os.path.join(build_path, 'dials', 'autocomplete', 'init.sh'))
    print

if __name__ == '__main__':
  gen = CompletionGenerator()
  gen.generate()
