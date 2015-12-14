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

from __future__ import division

try:
  # try importing scipy.linalg before any cctbx modules to avoid segfault on
  # some platforms
  import scipy.linalg # import dependency
except ImportError, e:
  pass

from cmd import Cmd
import sys

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


class ActionError(RuntimeError):
  '''
  Class to represent exception for when an action can't be performed

  '''

  def __init__(self, action, parent_action):
    '''
    :param action: The action we want to do
    :param parent_action: The parent's action

    '''
    text = 'Unable to perform "%s" after "%s"' % (action, parent_action)
    super(ActionError, self).__init__(text)


class ExternalCommand(object):
  '''
  Class to run an external command

  '''

  def __init__(self, command, filename=None, output=sys.stdout, wait_time=0.1):
    '''
    Run the command

    :param command: The command to run
    :param filename: The filename for stdout and stderr
    :param output: Write stdout and stderr to additional output
    :param wait_time: Wait time for polling file output

    '''
    import subprocess
    import time

    # Create the command string
    if not isinstance(command, str):
      command = subprocess.list2cmdline(command)


    # Open the file for output and input
    if filename is not None:
      with open(filename, "w") as outfile:
        with open(filename, "r") as infile:

          # Start the process
          process = subprocess.Popen(
            command,
            stdout=outfile,
            stderr=outfile,
            shell=True,
            universal_newlines=True,
            bufsize=-1)

          # Write the lines of the file to the output
          if output is not None:
            while True:
              line = infile.readline()
              if not line:
                if process.poll() is not None:
                  break
                time.sleep(wait_time)
              else:
                output.write(line)

          # Get the result
          self.result = process.wait()

    else:
      # Start the process
      process = subprocess.Popen(
        command,
        stdout=output,
        stderr=output,
        shell=True,
        universal_newlines=True,
        bufsize=-1)

      # Get the result
      self.result = process.wait()



def run_external_command(command, filename=None, output=sys.stdout, wait_time=0.1):
  '''
  Helper function to run command

  :param command: The command to run
  :param filename: The filename to output
  :param output: The buffer to write to
  :param wait_time: The polling timeout

  '''
  command = ExternalCommand(command, filename, output, wait_time)
  if command.result != 0:
    raise RuntimeError('Error: external command failed')


class ParameterManager(object):
  '''
  A class to manage the current set of parameters.

  '''

  def __init__(self, phil_scope):
    '''
    Create the master phil and set working phil to default parameters

    '''
    self.master_phil = phil_scope
    self.reset()

  def reset(self):
    '''
    Reset the working phil to the default parameters

    '''
    from libtbx.phil import parse
    self.working_phil = self.master_phil.fetch(source=parse(''))

  def set(self, parameters, short_syntax=False):
    '''
    Set a parameter and update the working phil
    :param parameter: The text string of parameters
    :param short_syntax: True/False treat as command line parameter

    '''
    from libtbx.phil import parse
    from libtbx.utils import Sorry
    import shlex
    if short_syntax == True:
      for parameter in shlex.split(parameters):
        interpretor = self.master_phil.command_line_argument_interpreter()
        self.working_phil = self.working_phil.fetch(
          interpretor.process_arg(parameter))
    else:
      working_phil, unused = self.working_phil.fetch(
        source=parse(parameters),
        track_unused_definitions=True)
      if len(unused) > 0:
        msg = [item.object.as_str().strip() for item in unused]
        msg = '\n'.join(['  %s' % line for line in msg])
        raise Sorry('The following definitions were not recognised\n%s' % msg)
      self.working_phil = working_phil

  def get(self, diff=True):
    '''
    Get the phil parameters

    :param diff: Get the diff phil

    '''
    if diff == False:
      result = self.working_phil
    else:
      result = self.master_phil.fetch_diff(source=self.working_phil)
    return result


class ImportParameterManager(ParameterManager):
  '''
  Specialization for import parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      include scope dials.command_line.import.phil_scope
    ''', process_includes=True)
    super(ImportParameterManager, self).__init__(phil_scope)


class FindSpotsParameterManager(ParameterManager):
  '''
  Specialization for find spots parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        datablock = None
          .type = str
      }
      include scope dials.command_line.find_spots.phil_scope
    ''', process_includes=True)
    super(FindSpotsParameterManager, self).__init__(phil_scope)


class DiscoverBetterModelParameterManager(ParameterManager):
  '''
  Specialization for discover_better_experiment_model parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        datablock = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.discover_better_experimental_model.phil_scope
    ''', process_includes=True)
    super(DiscoverBetterModelParameterManager, self).__init__(phil_scope)


class IndexParameterManager(ParameterManager):
  '''
  Specialization for index parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        datablock = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.index.phil_scope
    ''', process_includes=True)
    super(IndexParameterManager, self).__init__(phil_scope)


class RefineBSParameterManager(ParameterManager):
  '''
  Specialization for refine_bravais_settings parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.refine_bravais_settings.phil_scope
    ''', process_includes=True)
    super(RefineBSParameterManager, self).__init__(phil_scope)


class ReIndexParameterManager(ParameterManager):
  '''
  Specialization for reindex parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      solution = None
        .type = int
      input {
        reflections = None
          .type = str
      }
      include scope dials.command_line.reindex.phil_scope
    ''', process_includes=True)
    super(ReIndexParameterManager, self).__init__(phil_scope)


class RefineParameterManager(ParameterManager):
  '''
  Specialization for refine parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.refine.phil_scope
    ''', process_includes=True)
    super(RefineParameterManager, self).__init__(phil_scope)


class IntegrateParameterManager(ParameterManager):
  '''
  Specialization for integrate parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.integrate.phil_scope
    ''', process_includes=True)
    super(IntegrateParameterManager, self).__init__(phil_scope)


class ExportParameterManager(ParameterManager):
  '''
  Specialization for export parameters

  '''

  def __init__(self):
    '''
    Import phil scope and set up

    '''
    from libtbx.phil import parse
    phil_scope = parse('''
      description=None
        .type = str
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.export.phil_scope
    ''', process_includes=True)
    super(ExportParameterManager, self).__init__(phil_scope)


class GlobalParameterManager(dict):
  '''
  Class to hold all parameter managers

  '''

  def __init__(self):
    '''
    Init everything

    '''
    super(GlobalParameterManager, self).__init__()
    self.update({
      'import'                             : ImportParameterManager(),
      'find_spots'                         : FindSpotsParameterManager(),
      'discover_better_experimental_model' : DiscoverBetterModelParameterManager(),
      'index'                              : IndexParameterManager(),
      'refine_bravais_settings'            : RefineBSParameterManager(),
      'reindex'                            : ReIndexParameterManager(),
      'refine'                             : RefineParameterManager(),
      'integrate'                          : IntegrateParameterManager(),
      'export'                             : ExportParameterManager(),
    })


class Counter(object):
  '''
  A counter class to update command indices

  '''

  def __init__(self):
    '''
    Counter begins at zero

    '''
    self.count = 0

  def current(self):
    '''
    :return: The current counter

    '''
    return self.count

  def next(self):
    '''
    Update the counter value

    :return: The new counter value

    '''
    result = self.count
    self.count += 1
    return result


class CommandNode(object):
  '''
  A class to represent the commands

  '''

  parent_actions = []

  def __init__(self, parent=None, action='', parameters=None, directory=None):
    '''
    Initialise the tree parent and children and set the index

    :param parent: The command parent

    '''
    from os.path import join
    import copy

    # Check the parent is OK
    if parent is None:
      parent_action = None
    else:
      parent_action = parent.action
    if parent_action not in self.parent_actions:
      raise ActionError(action, parent_action)

    # Raise exception if trying to job after failure
    if parent is not None and parent.success == False:
      raise RuntimeError('Error: parent job %d failed' % parent.index)

    # Set the parent and counter
    self.parent = parent
    if self.parent is None:
      self.counter = Counter()
    else:
      self.counter = self.parent.counter

    # Set the index and some tree stuff
    self.index = self.counter.next()
    self.children = []
    if self.parent is not None:
      self.parent.children.append(self)

    # Initialise the description
    self.description = None

    # Init the result
    self.success = False

    # Save the info
    self.action = action
    self.params = copy.deepcopy(parameters)

    # Init the important paths
    if directory is not None:
      self.directory = join(directory, "%d_%s" % (self.index, self.action))
      self.output = join(self.directory, "output.txt")
      self.parameters = join(self.directory, "parameters.phil")
    else:
      self.directory = None
      self.output = None
      self.parameters = None
    self.report = None
    self.summary = None
    self.datablock = None
    self.experiments = None
    self.reflections = None

  def __iter__(self):
    '''
    Iterate through the children and their children

    '''
    yield self, 0
    for child in self.children:
      for node, depth in child:
        yield node, depth+1

  def apply(self):
    '''
    Apply the command

    :return: True/False success or failure

    '''
    from os.path import exists, join
    from os import makedirs

    # Check the output path does not exist already
    if exists(self.directory):
      raise RuntimeError('Output directory %s already exists' % self.directory)

    # Get the description
    self.description = self.params.get(diff=False).extract().description
    self.params.set("description=None")

    # Initialise running the command
    self.initialize()

    # Make the directory to store output
    makedirs(self.directory)

    # Set the parameter filename and write to file
    with open(self.parameters, "w") as outfile:
      outfile.write(self.params.get(diff=True).as_str())
    outfile.close()

    # Run the command (override this method)
    self.run()

    # Grab the result (override this method)
    self.finalize()

    # Set success
    self.success = True

  def generate_report(self):
    '''
    Helper function to run dials.report

    :param experiments: path to experiments.json
    :param reflections: path to reflections.json
    :param html: path to output html file

    '''
    command = ['dials.report']
    if self.reflections is None:
      raise RuntimeError('No reflections file set')
    if self.report is None:
      raise RuntimeError('No report file set')
    if self.experiments is not None:
      command.append('input.experiments=%s' % self.experiments)
    command.append('input.reflections=%s' % self.reflections)
    command.append('output.html=%s' % self.report)
    run_external_command(command)

  def check_files_exist(self, filenames=None):
    '''
    Helper function to check filenames exist

    '''
    from os.path import exists
    def assert_exists(name):
      if name is not None and name is not 'None' and not exists(name):
        raise RuntimeError("File %s could not be found" % name)
    if filenames is not None:
      for name in filenames:
        assert_exists(name)
    assert_exists(self.directory)
    assert_exists(self.output)
    assert_exists(self.parameters)
    assert_exists(self.report)
    assert_exists(self.summary)
    assert_exists(self.datablock)
    assert_exists(self.experiments)
    assert_exists(self.reflections)


class CommandTree(object):
  '''
  A class to provide to helpful tree functions

  '''

  def __init__(self, root):
    '''
    :param root: The tree root

    '''
    self.root = root

  def goto(self, index):
    '''
    Go to the desired node in the tree

    :param index: the index of the node to go to
    :return: The node at index

    '''
    for node, level in self.iternodes():
      if node.index == index:
        return node
    raise IndexError('Node %d not found' % index)

  def iternodes(self):
    '''
    Iterate through the tree nodes depth first

    '''
    for node, depth in self.root:
      yield node, depth

  def string(self, current=None):
    '''
    :return: The tree as a string

    '''
    size = len(str(self.root.counter.current()))
    def draw_tree(node, prefix):
      from cStringIO import StringIO
      buf = StringIO()
      if prefix:
        buf.write(('%%%dd' % size) % node.index)
        buf.write(' %s' % ('S' if node.success else 'F'))
        buf.write(prefix[:-3])
        buf.write('  +--')
      buf.write(node.action)
      if node.description is not None:
        buf.write(" %s" % node.description)
      if current is not None and node.index == current:
        buf.write(" (current)")
      buf.write('\n')
      for index, child in enumerate(node.children):
        if index+1 == len(node.children):
          sub_prefix = prefix + '   '
        else:
          sub_prefix = prefix + '  |'
        buf.write(draw_tree(child, sub_prefix))
      return buf.getvalue()
    return draw_tree(self.root, '')


class InitialState(CommandNode):
  '''
  A class to represent the initial clean state

  '''

  parent_actions = [None]

  def __init__(self):
    '''
    Initialise the command

    '''
    super(InitialState, self).__init__(None, 'clean')

    # Set success to True
    self.success = True

  def apply(self):
    '''
    Override apply since we have nothing to do here

    '''
    raise RuntimeError("Programming error: nothing to do")


class ImportCommand(CommandNode):
  '''
  A command to perform an import operation

  '''

  parent_actions = ['clean']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ImportCommand, self).__init__(
      parent, 'import', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.datablock = join(self.directory, "datablock.json")

    # Set filenames and input
    self.filenames = {
      'output.datablock' : self.datablock,
      'output.log'       : join(self.directory, "info.log"),
      'output.debug_log' : join(self.directory, "debug.log")
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the import command

    '''
    run_external_command(['dials.import', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Check the files exist
    self.check_files_exist(self.filenames.values())


class FindSpotsCommand(CommandNode):
  '''
  A command to perform an find_spots operation

  '''

  parent_actions = ['import']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(FindSpotsCommand, self).__init__(
      parent, 'find_spots', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.datablock = join(self.directory, "datablock.json")
    self.reflections = join(self.directory, "reflections.pickle")
    self.report = join(self.directory, 'report.html')

    # Set filenames and input
    self.filenames = {
      'input.datablock'    : self.parent.datablock,
      'output.datablock'   : self.datablock,
      'output.reflections' : self.reflections,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the find_spots command

    '''
    # Run find spots
    run_external_command(['dials.find_spots', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class DiscoverBetterModelCommand(CommandNode):
  '''
  A command to perform a discover_better_experimental_model operation

  '''

  parent_actions = ['find_spots']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(DiscoverBetterModelCommand, self).__init__(
      parent, 'discover_better_experimental_model', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.datablock = join(self.directory, "datablock.json")
    self.reflections = self.parent.reflections
    self.report = join(self.directory, 'report.html')

    # Set filenames and input
    self.filenames = {
      'input.datablock'    : self.parent.datablock,
      'input.reflections'  : self.parent.reflections,
      'output.datablock'   : self.datablock,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the index command

    '''
    run_external_command(['dials.discover_better_experimental_model',
                          self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class IndexCommand(CommandNode):
  '''
  A command to perform an index operation

  '''

  parent_actions = ['find_spots', 'discover_better_experimental_model']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(IndexCommand, self).__init__(
      parent, 'index', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.experiments = join(self.directory, "experiments.json")
    self.reflections = join(self.directory, "reflections.pickle")
    self.report = join(self.directory, 'report.html')

    # Set filenames and input
    self.filenames = {
      'input.datablock'    : self.parent.datablock,
      'input.reflections'  : self.parent.reflections,
      'output.reflections' : self.reflections,
      'output.experiments' : self.experiments,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the index command

    '''
    run_external_command(['dials.index', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class RefineBSCommand(CommandNode):
  '''
  A command to perform an refine_bravais_settings operation

  '''

  parent_actions = ['index']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(RefineBSCommand, self).__init__(
      parent, 'refine_bravais_settings', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.reflections = self.parent.reflections

    # Set other filenames
    self.filenames = {
      'input.experiments'  : self.parent.experiments,
      'input.reflections'  : self.parent.reflections,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.directory'   : self.directory,
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the refine_bravais_settings command

    '''
    run_external_command(['dials.refine_bravais_settings', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists, join
    import json

    # Read the summary and check all json files exist
    with open(join(self.directory, 'bravais_summary.json')) as summary_file:
      self.bs_summary = json.load(summary_file)
      self.bs_filenames = {}
      for name, value in self.bs_summary.iteritems():
        self.bs_filenames[name] = join(self.directory, 'bravais_setting_%s.json' % name)
      for name, value in self.bs_filenames.iteritems():
        if not exists(value):
          raise RuntimeError("File %s could not be found" % value)

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class ReIndexCommand(CommandNode):
  '''
  A command to perform an reindex operation

  '''

  parent_actions = ['refine_bravais_settings']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ReIndexCommand, self).__init__(
      parent, 'reindex', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # Get the solution we want and convert to the change_of_basis_op
    solution = self.params.get(diff=False).extract().solution
    if solution is None or solution == 'None':
      solution = max(map(int, self.parent.bs_summary.keys()))
    change_of_basis_op = self.parent.bs_summary[str(solution)]['cb_op']

    # Set the output experiments to the bravais settings file
    self.experiments = self.parent.bs_filenames[str(solution)]
    self.reflections = join(self.directory, "reflections.pickle")
    self.report = join(self.directory, 'report.html')

    # The files which can be set as parameters
    self.filenames = {
      'input.reflections'  : self.parent.reflections,
      'output.reflections' : self.reflections,
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

    # Set the solution parameter to None and set the cb_op
    self.params.set("solution=None")
    self.params.set("change_of_basis_op=%s" % change_of_basis_op)

  def run(self):
    '''
    Run the index command

    '''
    run_external_command(['dials.reindex', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class RefineCommand(CommandNode):
  '''
  A command to perform an refine operation

  '''

  parent_actions = ['index', 'reindex', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(RefineCommand, self).__init__(
      parent, 'refine', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    # set the results
    self.experiments = join(self.directory, "experiments.json")
    self.reflections = join(self.directory, "reflections.pickle")
    self.report = join(self.directory, 'report.html')

    # Set the other filenames and input
    self.filenames = {
      'input.experiments'  : self.parent.experiments,
      'input.reflections'  : self.parent.reflections,
      'output.reflections' : self.reflections,
      'output.experiments' : self.experiments,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.matches'     : join(self.directory, "matches.pickle"),
      'output.centroids'   : join(self.directory, "centroids.txt"),
      'output.history'     : join(self.directory, "history.txt"),
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the refine command

    '''
    run_external_command(['dials.refine', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class IntegrateCommand(CommandNode):
  '''
  A command to perform an integrate operation

  '''

  parent_actions = ['index', 'reindex', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(IntegrateCommand, self).__init__(
      parent, 'integrate', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join

    # set the results
    self.experiments = join(self.directory, "experiments.json")
    self.reflections = join(self.directory, "reflections.pickle")
    self.report = join(self.directory, 'report.html')

    # Set the other filenames and input
    self.filenames = {
      'input.experiments'  : self.parent.experiments,
      'input.reflections'  : self.parent.reflections,
      'output.reflections' : self.reflections,
      'output.experiments' : self.experiments,
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
      'output.report'      : join(self.directory, "summary.json"),
      'output.phil'        : 'None'
    }
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the integrate command

    '''
    run_external_command(['dials.integrate', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class ExportCommand(CommandNode):
  '''
  A command to perform an export operation

  '''

  parent_actions = ['integrate']

  def __init__(self, parent, parameters, directory):
    '''
    Initialise the command node

    :param parent: The parent command
    :param parameters: The parameters to use
    :param directory: The output directory

    '''
    super(ExportCommand, self).__init__(
      parent, 'export', parameters, directory)

  def initialize(self):
    '''
    Initialise the processing

    '''
    from os.path import join
    self.filenames = {
      'input.experiments'  : self.parent.experiments,
      'input.reflections'  : self.parent.reflections,
      'mtz.hklout'         : join(self.directory, "reflections.mtz"),
      'output.log'         : join(self.directory, "info.log"),
      'output.debug_log'   : join(self.directory, "debug.log"),
    }
    self.params.set("format=mtz")
    for name, value in self.filenames.iteritems():
      self.params.set('%s=%s' % (name, value))

  def run(self):
    '''
    Run the export command

    '''
    run_external_command(['dials.export', self.parameters], self.output)

  def finalize(self):
    '''
    Finalize the processing

    '''
    from os.path import exists
    import shutil

    # Check the files exist
    self.check_files_exist(self.filenames.values())

    # Copy the resulting mtz file to the working directory
    result_filename = "%d_integrated.mtz" % self.index
    shutil.copy2(self.filenames['mtz.hklout'], result_filename)


class ApplicationState(object):
  '''
  A class to hold all the application state

  '''

  def __init__(self, directory):
    '''
    Initialise the state

    :param directory: The output directory

    '''
    # Create the parameters
    self.parameters = GlobalParameterManager()

    # Set the initial state to current
    self.current = InitialState()

    # Create the command tree
    self.command_tree = CommandTree(self.current)

    # Save the parameters and directory
    self.directory = directory

    # Initialise the mode
    self.mode = 'import'

  def run(self):
    '''
    Run the command for the given mode

    '''

    # The command classes
    CommandClass = {
      'import'                             : ImportCommand,
      'find_spots'                         : FindSpotsCommand,
      'discover_better_experimental_model' : DiscoverBetterModelCommand,
      'index'                              : IndexCommand,
      'refine_bravais_settings'            : RefineBSCommand,
      'reindex'                            : ReIndexCommand,
      'refine'                             : RefineCommand,
      'integrate'                          : IntegrateCommand,
      'export'                             : ExportCommand
    }

    # Create the command
    command = CommandClass[self.mode](
      self.current,
      self.parameters[self.mode],
      self.directory)

    # Apply the command
    command.apply()

    # If successful update current
    self.current = command

  def goto(self, index):
    '''
    Goto a specific command

    :param index: The command index

    '''
    self.current = self.command_tree.goto(index)

  def history(self):
    '''
    Get the command history

    :return The command history

    '''
    return self.command_tree.string(current=self.current.index)

  def dump(self, filename):
    '''
    Dump the state to file

    :param filename: The filename

    '''
    import cPickle as pickle
    with open(filename, "w") as outfile:
      pickle.dump(self, outfile)

  @classmethod
  def load(Class, filename):
    '''
    Load the state from file

    :param filename: The filename
    :return: The state object

    '''
    import cPickle as pickle
    with open(filename) as infile:
      return pickle.load(infile)


class Controller(object):
  '''
  The controller class.

  This defines the interface the DIALS GUI and CLI programs can use to interact
  with the DIALS programs in a standard way.

  '''

  # The list of program modes
  mode_list = [
    'import',
    'find_spots',
    'discover_better_experimental_model',
    'index',
    'refine_bravais_settings',
    'reindex',
    'refine',
    'integrate',
    'export']

  def __init__(self,
               directory=".",
               state_filename="dials.state",
               recover=True):
    '''
    Initialise the controller

    :param directory: The output directory
    :param state_filename: The filename to save the state to
    :param recover: Recover the state if available

    '''
    from os.path import exists, abspath, join

    # Set some stuff
    self.state_filename = join(directory, state_filename)

    # Read state if available
    if recover == True and exists(state_filename):
      self.state = ApplicationState.load(state_filename)
      print "Recovered state from %s" % state_filename
      print self.get_history()
    else:
      def find_directory(working_directory):
        counter = 1
        while True:
          directory = join(working_directory, "dials-%d" % counter)
          if not exists(directory):
            return directory
          counter += 1
      self.state = ApplicationState(find_directory(abspath(directory)))

  def set_mode(self, mode):
    '''
    Set the current mode.

    :param mode: The mode to set

    '''
    # Is mode available?
    if mode not in self.mode_list:
      raise RuntimeError('Unknown mode: %s' % mode)

    # Set the mode
    self.state.mode = mode
    self.state.dump(self.state_filename)

  def get_mode(self):
    '''
    Get the current mode

    :return: The current mode

    '''
    return self.state.mode

  def set_parameters(self, parameters, short_syntax=False):
    '''
    Set the parameters.

    :param parameters: The parameters as a string
    :param show_syntax: Use command line string

    '''
    from libtbx.utils import Sorry
    self.state.parameters[self.get_mode()].set(parameters, short_syntax=short_syntax)
    self.state.dump(self.state_filename)

  def reset_parameters(self):
    '''
    Reset the parameters to the default values

    '''
    self.state.parameters[self.get_mode()].reset()
    self.state.dump(self.state_filename)

  def get_parameters(self, diff=True, mode=None):
    '''
    Get the current parameters

    :param diff: Show only the modified parameters

    '''
    if mode is None:
      mode = self.get_mode()
    return self.state.parameters[mode].get(diff=diff)

  def get_history(self):
    '''
    Get the history as a string

    :return: The history string

    '''
    return self.state.history()

  def get_models(self):
    '''
    Get the models filename

    :return: The models filename

    '''
    return self.state.current.models

  def get_report(self):
    '''
    Get the results filename

    :return: The results filename

    '''
    return self.state.current.report

  def get_summary(self):
    '''
    Get the report filename

    :return: The report filename

    '''
    return self.state.current.summary

  def goto(self, index):
    '''
    Change state to a different index

    :param index: The index to go to

    '''
    self.state.goto(index)
    self.state.dump(self.state_filename)

  def run(self):
    '''
    Run a program

    '''
    try:
      self.state.run()
      self.state.dump(self.state_filename)
    except Exception:
      self.state.dump(self.state_filename)
      raise


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

    # Initialise the console base
    Cmd.__init__(self)

    # Create the controller object
    self.controller = Controller()

    # Set the prompt to show the current mode
    self.prompt = "%s >> " % self.controller.get_mode()

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
    except Exception, e:
      print_error(e)

  def do_models(self, line):
    ''' Show the models '''
    import subprocess
    try:
      filename = self.controller.get_models()
      if filename is None:
        raise RuntimeError('No models to show')
      subprocess.call('dials.show %s' % filename, shell=True)
    except Exception, e:
      print_error(e)

  def do_summary(self, line):
    ''' Get the report. '''
    try:
      filename = self.controller.get_summary()
      if filename is None:
        raise RuntimeError('No result to show')
      print 'For report, see: %s' % filename
    except Exception, e:
      print_error(e)

  def do_report(self, line):
    ''' Get the results. '''
    try:
      import webbrowser
      filename = self.controller.get_report()
      if filename is None:
        raise RuntimeError('No result to show')
      webbrowser.open('file://%s' % filename)
    except Exception, e:
      print_error(e)

  def do_set(self, parameter):
    ''' Set a phil parameter '''
    try:
      self.controller.set_parameters(parameter, short_syntax=True)
    except Exception, e:
      print_error(e)

  def do_reset(self, line):
    ''' Reset parameters to default. '''
    try:
      self.controller.reset_parameters()
    except Exception, e:
      print_error(e)

  def do_load(self, filename):
    ''' Load a phil parameter file '''
    try:
      with open(filename) as infile:
        self.controller.set_parameters(infile.read())
    except Exception, e:
      print_error(e)

  def do_run(self, line):
    ''' Run a program '''
    try:
      self.controller.run()
      self.print_history()
    except Exception, e:
      print_error(e)

  def do_goto(self, line):
    ''' Goto a particular history state '''
    try:
      self.controller.goto(int(line))
      self.print_history()
    except Exception, e:
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
    except Exception, e:
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
      self.controller.run()
      self.print_history()
    except Exception, e:
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
