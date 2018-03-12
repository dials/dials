from __future__ import absolute_import, division, print_function

import os

try:
  from glob import glob
  import libtbx.load_env
  dials_path = libtbx.env.dist_path('dials')
  filenames = glob(os.path.join(dials_path, "extensions", "*.pyc"))
  if filenames:
    print("Cleaning up 'dials/extensions':")
    for filename in filenames:
      print(" Deleting %s" % filename)
      os.remove(filename)
except Exception:
  pass

try:
  from dials.framework import env
  env.cache.wipe()
except Exception:
  pass

try:
  from dials.util.version import dials_version
  print(dials_version())
except Exception:
  pass

import libtbx.pkg_utils
libtbx.pkg_utils.require('mock', '>=2.0')
libtbx.pkg_utils.require('orderedset')
libtbx.pkg_utils.require('pytest', '>=3.1')
libtbx.pkg_utils.require('Jinja2')
libtbx.pkg_utils.require('procrunner')

def _install_dials_autocompletion():
  '''generate bash.sh and SConscript file in /build/dials/autocomplete'''
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
  print('Identifying autocompletable commands:', end=' ')
  for file in sorted(os.listdir(commands_dir)):
    if not file.startswith('_') and file.endswith('.py'):
      if 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in open(os.path.join(commands_dir, file)).read():
        command_name = 'dials.%s' % file[:-3]
        print(command_name, end=' ')
        command_list.append(command_name)
  print()

  # Generate the autocompletion SConscript.
  with open(os.path.join(output_directory, 'SConscript'), 'w') as builder:
    builder.write('''Import("env")
import os.path
import libtbx.load_env
def dispatcher_outer(name):
  return os.path.join(libtbx.env.under_build('bin'), name)
def dispatcher_inner(name):
  return os.path.join(libtbx.env.dist_path('dials'), 'command_line', '%%s.py' %% name.partition('.')[2])
env.Append( BUILDERS={'AutoComplete': Builder(action='-$SOURCE --export-autocomplete-hints > $TARGET')} )
for cmd in [%s]:
  ac = env.AutoComplete(cmd, [dispatcher_outer(cmd), dispatcher_inner(cmd)])
  Requires(ac, Dir(libtbx.env.under_build('lib')))
  Depends(ac, os.path.join(libtbx.env.dist_path('dials'), 'util', 'options.py'))
  Depends(ac, os.path.join(libtbx.env.dist_path('dials'), 'util', 'autocomplete.sh'))
''' % ', '.join(["'%s'" % cmd for cmd in command_list]))

  # Generate a bash script activating command line completion for each relevant command
  with open(os.path.join(output_directory, 'bash.sh'), 'w') as script:
    script.write("type compopt &>/dev/null && {\n")
    for cmd in command_list:
      script.write(" complete -F _dials_autocomplete %s\n" % cmd)
    script.write("}\n")
    script.write("type compopt &>/dev/null || {\n")
    for cmd in command_list:
      script.write(" complete -o nospace -F _dials_autocomplete %s\n" % cmd)
    script.write("}\n")

  # Find the dials build directory
  build_path = abs(libtbx.env.build_path)

  # Permanently install the autocompletion script into setpaths-scripts.
  print("Installing autocompletion script into:", end=' ')
  for file in os.listdir(build_path):
    if file.startswith('setpath') and file.endswith('.sh'):
      original_file = open(os.path.join(build_path, file)).read()
      if not 'DIALS_ENABLE_COMMAND_LINE_COMPLETION' in original_file:
        marker = "\nexport PATH\n"
        original_position = original_file.find(marker)
        if original_position >= 0:
          print(file, end=' ')
          insert_position = original_position + len(marker)
          added_script = \
            '# DIALS_ENABLE_COMMAND_LINE_COMPLETION\n' \
            '[ -n "$BASH_VERSION" ] && {\n' \
            ' source $(libtbx.find_in_repositories dials/util/autocomplete.sh) && source %s || echo dials command line completion not available\n' \
            '}\n' % (
              os.path.join('$LIBTBX_BUILD', 'dials', 'autocomplete', 'bash.sh'))
          with open(os.path.join(build_path, file), 'w') as script:
            script.write(original_file[:insert_position] +
                         added_script +
                         original_file[insert_position:])
  print()

_install_dials_autocompletion()
