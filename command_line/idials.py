#!/usr/bin/env python
#
# dials.idials.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME idials

from __future__ import absolute_import, division

try:
  # try importing scipy.linalg before any cctbx modules to avoid segfault on
  # some platforms
  import scipy.linalg # import dependency
except ImportError:
  pass

from cmd import Cmd

help_message = '''

This program is an interactive console mode for the dials command line programs.
It provides convient access to the following dials programs

 - dials.import
 - dials.find_spots
 - dials.discover_better_experimental_model
 - dials.index
 - dials.refine_bravais_settings
 - dials.reindex
 - dials.refine
 - dials.integrate
 - dials.export

To view the full list of commands available in the console, type "help".

The program implements a tree model for the data processing. Each command
executed is treated as a node in the tree and each job is output into its own
directory where all output associated with the job will be kept. The current
state can be changed by using the "goto" command to jump between nodes of the
tree. This allows the user to go back and re-do a stage of the processing if
necessary.

The program uses a modal state model. To change mode use the "mode" command.
Parameters are set using the "set" command. Parameters can be viewed using the
"get" command which shows only the modified parameters and the "all" command
which shows all the parameters. Parameters can be set to default values by using
the "reset" command. The mode is then run using the "run" command. This design
allows complicated parameter definitions to be specified in a more convenient
way for the user.

For convenience when the user wants to use default or a small number of
parameters, an imperative mode is supplied for the dials programs. Each of these
commands changes the mode, sets the parameters and runs the program in a single
line.

The program implements a persistant state mechanism. If the program crashes, the
state can be recovered by simply restarting the program in the same directory.

Examples::

  # Process data, select bravais setting and integrate with/without profile
  # fitting
  >> import template=/path/to/images_####.cbf
  >> find_spots
  >> index
  >> refine_bravais_settings
  >> reindex solution=1
  >> refine scan_varying=True
  >> integrate
  >> export
  >> goto 6
  >> integrate profile.fitting=False
  >> export

  # Simple scripting can be done by inputting command on stdin
  idials <<EOF
  import template=/path/to/images_####.cbf
  find_spots
  index
  refine scan_varying=True
  integrate
  export
  EOF

'''

def print_error(exception):
  '''
  Print out the error message

  '''
  print ''
  print '*' * 80
  print 'USER ERROR: PLEASE REPLACE USER'
  print ''
  print exception
  print '*' * 80
  print ''


class Console(Cmd):
  '''
  A class to implement an interactive dials console

  '''

  # The default prompt
  prompt = ">> "

  def __init__(self):
    '''
    Initialise the console

    '''
    from dials.util.idials import Controller

    # Initialise the console base
    Cmd.__init__(self)

    # Open a file for history
    self.command_history = open("command_history", "a")

    # Create the controller object
    self.controller = Controller()

    # Set the prompt to show the current mode
    self.prompt = "%s >> " % self.controller.get_mode()

  def precmd(self, line):
    '''
    Process command before

    '''
    # Add command history to file
    self.command_history.write("%s\n" % line)
    self.command_history.flush()

    # Call parent
    return Cmd.precmd(self, line)

  def emptyline(self):
    ''' Do nothing on empty line '''
    pass

  def default(self, line):
    ''' The default line handler '''
    try:
      self.controller.set_parameters(line, short_syntax=True)
    except Exception:
      return Cmd.default(self, line)

  def do_mode(self, mode):
    ''' Set the program mode '''
    try:
      self.controller.set_mode(mode)
      self.prompt = "%s >> " % self.controller.get_mode()
    except Exception as e:
      print_error(e)

  def do_models(self, line):
    ''' Show the models '''
    import subprocess
    try:
      filename = self.controller.get_models()
      if filename is None:
        raise RuntimeError('No models to show')
      subprocess.call('dials.show %s' % filename, shell=True)
    except Exception as e:
      print_error(e)

  def do_summary(self, line):
    ''' Get the report. '''
    try:
      filename = self.controller.get_summary()
      if filename is None:
        raise RuntimeError('No result to show')
      print 'For report, see: %s' % filename
    except Exception as e:
      print_error(e)

  def do_report(self, line):
    ''' Get the results. '''
    try:
      import webbrowser
      filename = self.controller.get_report()
      if filename is None:
        raise RuntimeError('No result to show')
      webbrowser.open('file://%s' % filename)
    except Exception as e:
      print_error(e)

  def do_set(self, parameter):
    ''' Set a phil parameter '''
    try:
      self.controller.set_parameters(parameter, short_syntax=True)
    except Exception as e:
      print_error(e)

  def do_reset(self, line):
    ''' Reset parameters to default. '''
    try:
      self.controller.reset_parameters()
    except Exception as e:
      print_error(e)

  def do_undo(self, line):
    ''' Undo parameters. '''
    try:
      self.controller.undo_parameters()
      print self.controller.get_parameters(diff=True).as_str()
    except Exception as e:
      print_error(e)

  def do_redo(self, line):
    ''' Redo parameters. '''
    try:
      self.controller.redo_parameters()
      print self.controller.get_parameters(diff=True).as_str()
    except Exception as e:
      print_error(e)

  def do_load(self, filename):
    ''' Load a phil parameter file '''
    try:
      with open(filename) as infile:
        self.controller.set_parameters(infile.read())
    except Exception as e:
      print_error(e)

  def do_run(self, line):
    ''' Run a program '''
    try:
      self.controller.run().wait()
      self.print_history()
    except Exception as e:
      print_error(e)

  def do_goto(self, line):
    ''' Goto a particular history state '''
    try:
      self.controller.goto(int(line))
      self.print_history()
    except Exception as e:
      print_error(e)

  def do_get(self, line):
    ''' Show all the possible parameters '''
    print self.controller.get_parameters(diff=True).as_str()

  def do_all(self, line):
    ''' Show all the possible parameters '''
    print self.controller.get_parameters(diff=False).as_str()

  def do_history(self, line):
    ''' Show the history. '''
    self.print_history()

  def do_import(self, line):
    ''' Imperative import command '''
    self.run_import_as_imperative(line)

  def do_find_spots(self, params):
    ''' Imperative find_spots command '''
    self.run_as_imperative("find_spots", params)

  def do_discover_better_experimental_model(self, params):
    ''' Imperative discover_better_experimental_model command '''
    self.run_as_imperative("discover_better_experimental_model", params)

  def do_index(self, params):
    ''' Imperative index command '''
    self.run_as_imperative("index", params)

  def do_refine_bravais_settings(self, params):
    ''' Imperative refine_bravais_settings command '''
    self.run_as_imperative("refine_bravais_settings", params)

  def do_reindex(self, params):
    ''' Imperative reindex command '''
    self.run_as_imperative("reindex", params)

  def do_refine(self, params):
    ''' Imperative refine command '''
    self.run_as_imperative("refine", params)

  def do_integrate(self, params):
    ''' Imperative integrate command '''
    self.run_as_imperative("integrate", params)

  def do_export(self, params):
    ''' Imperative export command '''
    self.run_as_imperative("export", params)

  def do_exit(self, line):
    ''' Exit the console '''
    return True

  def do_EOF(self, line):
    ''' Exit the console '''
    print ''
    return True

  def do_shell(self, line):
    ''' Execute shell commands '''
    import subprocess
    try:
      subprocess.call(line, shell=True)
    except Exception as e:
      print_error(e)

  def run_import_as_imperative(self, line):
    '''
    Helper for import imperative mode. Change mode, set parameters and run the job

    '''
    self.run_as_imperative("import", self.parse_import_line(line))

  def run_as_imperative(self, mode, parameters):
    '''
    Helper for imperative mode. Change mode, set parameters and run the job

    '''
    try:
      self.controller.set_mode(mode)
      self.prompt = "%s >> " % self.controller.get_mode()
      self.controller.set_parameters(parameters, short_syntax=True)
      self.controller.run().wait()
      self.controller.undo_parameters()
      self.print_history()
    except Exception as e:
      print_error(e)

  def parse_import_line(self, line):
    ''' Given a line after the import command. Figure out phil and filenames

    Split line like a shell command line. Then check each argument. If an
    argument is a directory, then find templates in that directory. Otherwise,
    Find files matching the argument using the glob module and generate
    templates from the matches. If there are no matches then input as a phil
    parameter.

    '''
    from os import listdir
    from os.path import isdir, isfile, join
    from glob import glob
    from dxtbx.model.scan_helpers import template_regex
    import shlex
    def templates_from_filenames(filenames):
      templates = []
      for f in filenames:
        try:
          templates.append(template_regex(f)[0])
        except Exception:
          pass
      return list(set(templates))
    def templates_from_directory(directory):
      filenames = []
      for f in listdir(directory):
        if isfile(join(directory, f)):
          filenames.append(join(directory, f))
      return templates_from_filenames(sorted(filenames))
    parameters = []
    arguments = shlex.split(line)
    for arg in arguments:
      if isdir(arg):
        templates = templates_from_directory(arg)
        for t in templates:
          parameters.append("template=%s" % t)
      else:
        matches = glob(arg)
        if len(matches) > 0:
          templates = templates_from_filenames(matches)
          for t in templates:
            parameters.append("template=%s" % t)
        else:
          parameters.append(arg)
    return ' '.join(parameters)

  def print_history(self):
    '''
    Print the history

    '''
    print ''
    print 'History'
    print self.controller.get_history()

  def complete_mode(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for changing mode.

    '''
    return [i for i in self.controller.mode_list if i.startswith(text)]

  def complete_set(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for setting parameters

    '''
    return self.get_phil_completions(text)

  def complete_load(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for loading phil files

    '''
    return get_path_completions(text, line, begidx, endidx)

  def complete_import(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for import

    '''
    return self.get_phil_completions(text, mode="import")


  def complete_find_spots(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for find_spots

    '''
    return self.get_phil_completions(text, mode="find_spots")

  def complete_discover_better_experimental_model(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for discover_better_experimental_model

    '''
    return self.get_phil_completions(text, mode="discover_better_experimental_model")

  def complete_index(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for index

    '''
    return self.get_phil_completions(text, mode="index")

  def complete_refine_bravais_settings(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for refine_bravais_settings

    '''
    return self.get_phil_completions(text, mode="refine_bravais_settings")

  def complete_reindex(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for reindex

    '''
    return self.get_phil_completions(text, mode="reindex")

  def complete_refine(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for refine

    '''
    return self.get_phil_completions(text, mode="refine")


  def complete_integrate(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for integrate

    '''
    return self.get_phil_completions(text, mode="integrate")

  def complete_export(self, text, line, begidx, endidx):
    '''
    Offer tab completion options for export

    '''
    return self.get_phil_completions(text, mode="export")

  def get_phil_completions(self, text, mode=None):
    '''
    Get completions for phil parameters

    '''
    phil_scope = self.controller.get_parameters(diff=False, mode=mode)
    definitions = phil_scope.all_definitions()
    full_names = [d.path for d in definitions]
    return [i for i in full_names if text in i]

  def get_path_completions(self, text, line, begidx, endidx):
    '''
    Get completions for paths

    '''
    import os
    from os.path import isdir
    import glob

    def _append_slash_if_dir(p):
      if p and isdir(p) and p[-1] != os.sep:
        return p + os.sep
      else:
        return p

    before_arg = line.rfind(" ", 0, begidx)
    if before_arg == -1:
      return # arg not found

    fixed = line[before_arg+1:begidx]  # fixed portion of the arg
    arg = line[before_arg+1:endidx]
    pattern = arg + '*'

    completions = []
    for path in glob.glob(pattern):
      path = _append_slash_if_dir(path)
      completions.append(path.replace(fixed, "", 1))
    return completions


# The intre string for the console
CONSOLE_INTRO = '''
DIALS interactive mode
Type "help" for more information
'''


if __name__ == '__main__':

  # Print the console intro
  print CONSOLE_INTRO

  # Create the console
  console = Console()

  # Enter the command loop
  console.cmdloop()
