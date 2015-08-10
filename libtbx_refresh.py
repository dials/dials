from __future__ import division

try:
  from glob import glob
  import libtbx.load_env
  import os
  dials_path = libtbx.env.dist_path('dials')
  filenames = glob(os.path.join(dials_path, "extensions", "*.pyc"))
  if len(filenames) > 0:
    print "Cleaning up 'dials/extensions':"
    for filename in filenames:
      print " Deleting %s" % filename
      os.remove(filename)
except Exception:
  pass

try:
  from dials.framework import env
  env.cache.wipe()
except Exception:
  pass

def _prepare_dials_autocompletion():
  '''generate init.sh and SConscript file in /build/dials/autocomplete'''
  import libtbx.load_env
  import os

  # Find the dials source directory
  dist_path = libtbx.env.dist_path('dials')

  # Set the location of the output directory
  output_directory = libtbx.env.under_build(os.path.join('dials', 'autocomplete'))
  try:
    os.makedirs(output_directory)
  except OSError:
    for file in os.listdir(output_directory):
      file_path = os.path.join(output_directory, file)
      try:
        if os.path.isfile(file_path):
          os.unlink(file_path)
      except Exception, e:
        pass
    pass

  commands_dir = os.path.join(dist_path, 'command_line')
  command_list = []
  print 'Identifying autocompletable commands:',
  for file in sorted(os.listdir(commands_dir)):
    if not file.startswith('_') and file.endswith('.py'):
      if 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in open(os.path.join(commands_dir, file)).read():
        command_name = 'dials.%s' % file[:-3]
        print command_name,
        command_list.append(command_name)
  print

  # Generate the autocompletion SConscript.
  with open(os.path.join(output_directory, 'SConscript'), 'w') as builder:
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

  # Generate the autocompletion bash init script
  with open(os.path.join(output_directory, 'init.sh'), 'w') as loader:
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
      os.path.join(output_directory, '${cmd}'),
      os.path.join(output_directory, '${cmd}'),
      os.path.join(output_directory, '${cmd}'),
      dist_path,
      "\n".join(["complete -F _dials_autocomplete %s" % cmd for cmd in command_list])
   ))

def _install_dials_autocompletion():
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

_prepare_dials_autocompletion()
#_install_dials_autocompletion()
