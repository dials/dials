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
    builder.write('''Import("env")
import os.path
import libtbx.load_env
def dispatcher_outer(name):
  return os.path.join(libtbx.env.under_build('bin'), name)
def dispatcher_inner(name):
  return os.path.join(libtbx.env.dist_path('dials'), 'command_line', '%%s.py' %% name.partition('.')[2])
env.Append( BUILDERS={'AutoComplete': Builder(action='$SOURCE --export-autocomplete-hints > $TARGET')} )
for cmd in [%s]:
  ac = env.AutoComplete(cmd, [dispatcher_outer(cmd), dispatcher_inner(cmd)])
  Requires(ac, Dir(libtbx.env.under_build('lib')))
''' % ( ', '.join(["'%s'" % cmd for cmd in command_list]) ))

  return command_list

def _install_dials_autocompletion(command_list):
  '''Permanently install the autocompletion script into setpaths-scripts.'''
  import libtbx.load_env
  import os

  # Find the dials source directory
  dist_path = libtbx.env.dist_path('dials')

  # Find the dials build directory
  build_path = abs(libtbx.env.build_path)
  print "Installing autocompletion script into:",
  for file in os.listdir(build_path):
    if file.startswith('setpath') and file.endswith('.sh'):
      original_file = open(os.path.join(build_path, file)).read()
      if not 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in original_file:
        marker = "\nexport PATH\n"
        insert_position = original_file.find(marker) + len(marker)
        if insert_position >= 0:
          print file,
          added_script = \
            '# DIALS_ENABLE_COMMAND_LINE_COMPLETION\n' \
            '[ -z "$BASH_VERSIONINFO" ] && {\n' \
            ' source %s\n' \
            ' %s\n' \
            '}\n' % (
              os.path.join('$LIBTBX_BUILD', '..', 'modules', 'dials', 'util', 'autocomplete.sh'),
              "\n ".join(["complete -F _dials_autocomplete %s" % cmd for cmd in command_list]))
          with open(os.path.join(build_path, file), 'w') as script:
            script.write(original_file[:insert_position] +
                         added_script + 
                         original_file[insert_position:])
  print

_autocomplete_commands = _prepare_dials_autocompletion()
#_install_dials_autocompletion(_autocomplete_commands)
