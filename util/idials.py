#!/usr/bin/env python
#
# idials.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import sys

try:
  # try importing scipy.linalg before any cctbx modules to avoid segfault on
  # some platforms
  import scipy.linalg # import dependency
except ImportError, e:
  pass


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



class UndoStack(object):
  '''
  A class to implement an undo stack

  '''

  def __init__(self, initial):
    '''
    Set the initial state

    '''
    self._stack = [initial]
    self._index = 0

  def push(self, obj):
    '''
    Add a new state

    '''
    self._stack = self._stack[0:self._index+1] + [obj]
    self._index = len(self._stack) - 1

  def peek(self):
    '''
    Get the current state

    '''
    return self._stack[self._index]

  def undo(self):
    '''
    Undo to the last state

    '''
    if self._index > 0:
      self._index -= 1

  def redo(self):
    '''
    Redo to the next state

    '''
    if self._index < len(self._stack) - 1:
      self._index += 1

  def reset(self):
    '''
    Reset to the initial state

    '''
    self._stack = self._stack[0:1]
    self._index = 0

  def __len__(self):
    '''
    Get the number of states

    '''
    return len(self._stack)


class ParameterManager(object):
  '''
  A class to manage the current set of parameters.

  '''

  def __init__(self, phil_scope):
    '''
    Create the master phil and set working phil to default parameters

    '''
    from libtbx.phil import parse
    self.master_phil = phil_scope
    self.working_phil = UndoStack(self.master_phil.fetch(source=parse('')))

  def undo(self):
    '''
    Undo the last parameter changes

    '''
    self.working_phil.undo()

  def redo(self):
    '''
    Redo the last parameter changes

    '''
    self.working_phil.redo()

  def reset(self):
    '''
    Reset the working phil to the default parameters

    '''
    from libtbx.phil import parse
    self.working_phil.push(self.master_phil.fetch(source=parse('')))

  def set(self, parameters, short_syntax=False):
    '''
    Set a parameter and update the working phil
    :param parameter: The text string of parameters
    :param short_syntax: True/False treat as command line parameter

    '''
    from libtbx.phil import parse
    from libtbx.utils import Sorry
    import shlex
    if short_syntax:
      working_phil = self.working_phil.peek()
      for parameter in shlex.split(parameters):
        interpretor = self.master_phil.command_line_argument_interpreter()
        working_phil = working_phil.fetch(
          interpretor.process_arg(parameter))
    else:
      working_phil, unused = self.working_phil.peek().fetch(
        source=parse(parameters),
        track_unused_definitions=True)
      if len(unused) > 0:
        msg = [item.object.as_str().strip() for item in unused]
        msg = '\n'.join(['  %s' % line for line in msg])
        raise Sorry('The following definitions were not recognised\n%s' % msg)
    self.working_phil.push(working_phil)

  def get(self, diff=True):
    '''
    Get the phil parameters

    :param diff: Get the diff phil

    '''
    working_phil = self.working_phil.peek()
    if diff:
      result = self.master_phil.fetch_diff(source=working_phil)
    else:
      result = working_phil
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
        experiments = None
          .type = str
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
    phil_scope = phil_scope.fetch(parse("mtz.hklout=None"))
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

  def __init__(self, count=1):
    '''
    Counter begins at zero

    '''
    self.count = count

  def current(self):
    '''
    :return: The current counter

    '''
    return self.count

  def incr(self):
    '''
    Update the counter value

    '''
    self.count += 1


class CommandState(object):
  '''
  A class to represent the command state

  '''
  def __init__(self, **kwargs):
    '''
    Initialise the command state

    '''
    self.applied = False
    self.success = False
    self.parent = None
    self.children = []
    self.__dict__.update(kwargs)

  def __iter__(self):
    '''
    Iterate through the children and their children

    '''
    yield self, 0
    for child in self.children:
      for node, depth in child:
        yield node, depth+1

  def as_dict(self):
    '''
    Return the command state as a dictionary

    '''
    dictionary = {}
    for key, value in self.__dict__.iteritems():
      if key == 'parent':
        continue
      elif key == 'children':
        value = [c.as_dict() for c in value]
      dictionary[key] = value
    return dictionary

  @classmethod
  def from_dict(cls, dictionary):
    '''
    Convert the dictionary to the state

    '''
    state = cls()
    for key, value in dictionary.iteritems():
      if key == 'children':
        for item in value:
          child = cls.from_dict(item)
          child.parent = state
          state.children.append(child)
      else:
        setattr(state, key, value)
    return state


class Command(object):
  '''
  A class to represent the commands

  '''

  name = None

  def __init__(self, parent=None, index=None, phil_scope=None, workspace=None):
    '''
    Initialise the action

    '''
    from os.path import join
    import copy

    # Check some input
    if self.name is None:
      raise RuntimeError('No action name set')
    if index is None:
      raise RuntimeError('No index set')
    if phil_scope is None:
      raise RuntimeError('No parameter set')
    if workspace is None:
      raise RuntimeError('No workspace set')

    # Check the parent
    if parent is not None:
      if not parent.applied:
        raise RuntimeError('Parent job %d not applied' % parent.index)
      if not parent.success:
        raise RuntimeError('Parent job %d failed' % parent.index)
      if parent.index >= index:
        raise RuntimeError('Invalid parent/child indices: %d / %d' % (parent.index, index))
      if parent.workspace is not None and parent.workspace != workspace:
        raise RuntimeError('Invalid parent/child worksapce: %s / %s' % (parent.workspace, workspace))
      if parent.name not in self.allowed_parents:
        raise RuntimeError('Invalid parent/child: %s / %s' % (parent.name, self.name))

    # Set the phil stuff and description
    description = phil_scope.get(diff=False).extract().description
    self.phil_scope = copy.deepcopy(phil_scope)
    self.phil_scope.set("description=None")

    # The directory
    directory = "%d_%s" % (index, self.name)

    # Set the state
    self.state = CommandState(
      name        = self.name,
      index       = index,
      parent      = parent,
      description = description,
      workspace   = workspace,
      directory   = join(workspace, directory),
      output      = join(workspace, directory, "output.txt"),
      parameters  = join(workspace, directory, "parameters.phil"))

    # Set the parent and children
    if self.state.parent is not None:
      self.state.parent.children.append(self.state)

  def apply(self):
    '''
    Apply the command

    '''
    from os.path import exists, join
    from os import makedirs

    # Check that the command has not already been run
    if self.state.applied:
      raise RuntimeError('This command has already been run')
    self.state.applied = True

    # Check the output path does not exist already
    if exists(self.state.directory):
      raise RuntimeError('Output directory %s already exists' % self.state.directory)

    # Make the directory to store output
    makedirs(self.state.directory)

    # Run the command (override this method)
    self.run()

    # Set success
    self.state.success = True

    # Return the state
    return self.state

  def run(self):
    '''
    Run the command

    '''
    pass

  def generate_report(self):
    '''
    Helper function to run dials.report

    :param experiments: path to experiments.json
    :param reflections: path to reflections.json
    :param html: path to output html file

    '''
    command = ['dials.report']
    if not hasattr(self.state, "reflections") or self.state.reflections is None:
      raise RuntimeError('No reflections file set')
    if not hasattr(self.state, "report") or self.state.report is None:
      raise RuntimeError('No report file set')
    if hasattr(self.state, "experiments") and self.state.experiments is not None:
      command.append('input.experiments=%s' % self.state.experiments)
    command.append('input.reflections=%s' % self.state.reflections)
    command.append('output.html=%s' % self.state.report)
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
    def assert_exists_if_present(name):
      if hasattr(self.state, name):
        assert_exists(getattr(self.state, name))
    assert_exists_if_present("directory")
    assert_exists_if_present("output")
    assert_exists_if_present("parameters")
    assert_exists_if_present("report")
    assert_exists_if_present("summary")
    assert_exists_if_present("datablock")
    assert_exists_if_present("experiments")
    assert_exists_if_present("reflections")


class CommandTree(object):
  '''
  A class to provide to helpful tree functions

  '''

  def __init__(self, root, counter):
    '''
    :param root: The tree root

    '''
    self.root = root
    self.counter = counter

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
    size = len(str(self.counter.current()))
    def draw_tree(node, prefix):
      from cStringIO import StringIO
      buf = StringIO()
      if prefix:
        buf.write(('%%%dd' % size) % node.index)
        buf.write(' %s' % ('S' if node.success else 'F'))
        buf.write(prefix[:-3])
        buf.write('  +--')
      buf.write(node.name)
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


class InitialState(CommandState):
  '''
  The initial state

  '''

  def __init__(self):
    '''
    Initialise the state

    '''
    super(InitialState, self).__init__(
      name        = "clean",
      index       = 0,
      description = '',
      workspace   = None,
      applied     = True,
      success     = True)

  def apply(self):
    '''
    Override apply since we have nothing to do here

    '''
    raise RuntimeError("Programming error: nothing to do")


class Import(Command):
  '''
  A command to perform an import operation

  '''

  name = 'import'

  allowed_parents = ['clean']

  def run(self):
    '''
    Run the import command

    '''
    from os.path import join

    # set the results
    self.state.datablock = join(self.state.directory, "datablock.json")

    # Set filenames and input
    self.filenames = {
      'output.datablock' : self.state.datablock,
      'output.log'       : join(self.state.directory, "info.log"),
      'output.debug_log' : join(self.state.directory, "debug.log")
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.import', self.state.parameters], self.state.output)

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class FindSpots(Command):
  '''
  A command to perform an find_spots operation

  '''

  name='find_spots'

  allowed_parents = ['import']

  def run(self):
    '''
    Run the find_spots command

    '''
    from os.path import join

    # set the results
    self.state.datablock = join(self.state.directory, "datablock.json")
    self.state.reflections = join(self.state.directory, "reflections.pickle")
    self.state.report = join(self.state.directory, 'report.html')

    # Set filenames and input
    self.filenames = {
      'input.datablock'    : self.state.parent.datablock,
      'output.datablock'   : self.state.datablock,
      'output.reflections' : self.state.reflections,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run find spots
    run_external_command(['dials.find_spots', self.state.parameters],
                         self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class DiscoverBetterExperimentalModel(Command):
  '''
  A command to perform a discover_better_experimental_model operation

  '''

  name='discover_better_experimental_model'

  allowed_parents = ['find_spots']

  def run(self):
    '''
    Run the index command

    '''
    from os.path import join

    # set the results
    self.state.datablock = join(self.state.directory, "datablock.json")
    self.state.reflections = self.state.parent.reflections
    self.state.report = join(self.state.directory, 'report.html')

    # Set filenames and input
    self.filenames = {
      'input.datablock'    : self.state.parent.datablock,
      'input.reflections'  : self.state.parent.reflections,
      'output.datablock'   : self.state.datablock,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.discover_better_experimental_model',
                          self.state.parameters], self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class Index(Command):
  '''
  A command to perform an index operation

  '''

  name = 'index'

  allowed_parents = [
    'find_spots',
    'discover_better_experimental_model',
    'index'
  ]

  def run(self):
    '''
    Run the index command

    '''
    from os.path import join

    # set the results
    self.state.experiments = join(self.state.directory, "experiments.json")
    self.state.reflections = join(self.state.directory, "reflections.pickle")
    self.state.report = join(self.state.directory, 'report.html')

    # Set filenames and input
    if self.state.parent.name == 'index':
      self.filenames = {
        'input.experiments' : self.state.parent.experiments,
      }
    else:
      self.filenames = {
        'input.datablock' : self.state.parent.datablock,
      }
    self.filenames.update({
      'input.reflections'  : self.state.parent.reflections,
      'output.reflections' : self.state.reflections,
      'output.experiments' : self.state.experiments,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
    })
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.index', self.state.parameters],
                         self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class RefineBravaisSettings(Command):
  '''
  A command to perform an refine_bravais_settings operation

  '''

  name='refine_bravais_settings'

  allowed_parents = ['index']

  def run(self):
    '''
    Run the refine_bravais_settings command

    '''
    from os.path import exists
    from os.path import join

    # set the results
    self.state.reflections = self.state.parent.reflections
    self.state.summary = join(self.state.directory, "bravais_summary.json")

    # Set other filenames
    self.filenames = {
      'input.experiments'  : self.state.parent.experiments,
      'input.reflections'  : self.state.parent.reflections,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
      'output.directory'   : self.state.directory,
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.refine_bravais_settings',
                          self.state.parameters], self.state.output)

    # Read the summary and check all json files exist
    for item, filename in self.bravais_setting_filenames(self.state).iteritems():
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)

    # Check the files exist
    self.check_files_exist(self.filenames.values())

  @classmethod
  def bravais_summary(cls, state):
    '''
    Get the bravais summary

    '''
    import json
    with open(state.summary) as summary_file:
      return json.load(summary_file)

  @classmethod
  def bravais_setting_filenames(cls, state):
    '''
    Get the bravais setting filenames

    '''
    from os.path import join
    bs_summary = cls.bravais_summary(state)
    bs_filenames = {}
    for name, value in bs_summary.iteritems():
      bs_filenames[name] = join(state.directory, 'bravais_setting_%s.json' % name)
    return bs_filenames


class Reindex(Command):
  '''
  A command to perform an reindex operation

  '''

  name='reindex'

  allowed_parents = ['refine_bravais_settings']

  def run(self):
    '''
    Run the index command

    '''
    from os.path import join

    # Get the solution we want and convert to the change_of_basis_op
    summary = RefineBravaisSettings.bravais_summary(self.state.parent)
    solution = self.phil_scope.get(diff=False).extract().solution
    if solution is None or solution == 'None':
      solution = max(map(int, summary.keys()))
    change_of_basis_op = summary[str(solution)]['cb_op']

    # Set the output experiments to the bravais settings file
    bs_filenames = RefineBravaisSettings.bravais_setting_filenames(self.state.parent)
    self.state.experiments = bs_filenames[str(solution)]
    self.state.reflections = join(self.state.directory, "reflections.pickle")
    self.state.report = join(self.state.directory, 'report.html')

    # The files which can be set as parameters
    self.filenames = {
      'input.reflections'  : self.state.parent.reflections,
      'output.reflections' : self.state.reflections,
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the solution parameter to None and set the cb_op
    self.phil_scope.set("solution=None")
    self.phil_scope.set("change_of_basis_op=%s" % change_of_basis_op)

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.reindex', self.state.parameters],
                         self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class Refine(Command):
  '''
  A command to perform an refine operation

  '''

  name = 'refine'

  allowed_parents = [
    'index',
    'reindex',
    'refine',
    'integrate'
  ]

  def run(self):
    '''
    Run the refine command

    '''
    from os.path import join

    # set the results
    self.state.experiments = join(self.state.directory, "experiments.json")
    self.state.reflections = join(self.state.directory, "reflections.pickle")
    self.state.report = join(self.state.directory, 'report.html')

    # Set the other filenames and input
    self.filenames = {
      'input.experiments'  : self.state.parent.experiments,
      'input.reflections'  : self.state.parent.reflections,
      'output.reflections' : self.state.reflections,
      'output.experiments' : self.state.experiments,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
      'output.matches'     : join(self.state.directory, "matches.pickle"),
      'output.centroids'   : join(self.state.directory, "centroids.pickle"),
      'output.history'     : join(self.state.directory, "history.pickle"),
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.refine', self.state.parameters],
                         self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class Integrate(Command):
  '''
  A command to perform an integrate operation

  '''

  name = 'integrate'

  allowed_parents = [
    'index',
    'reindex',
    'refine',
    'integrate'
  ]

  def run(self):
    '''
    Run the integrate command

    '''
    from os.path import join

    # set the results
    self.state.experiments = join(self.state.directory, "experiments.json")
    self.state.reflections = join(self.state.directory, "reflections.pickle")
    self.state.report = join(self.state.directory, 'report.html')

    # Set the other filenames and input
    self.filenames = {
      'input.experiments'  : self.state.parent.experiments,
      'input.reflections'  : self.state.parent.reflections,
      'output.reflections' : self.state.reflections,
      'output.experiments' : self.state.experiments,
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
      'output.report'      : join(self.state.directory, "summary.json"),
      'output.phil'        : 'None'
    }
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.integrate', self.state.parameters],
                         self.state.output)

    # Generate the report
    self.generate_report()

    # Check the files exist
    self.check_files_exist(self.filenames.values())


class Export(Command):
  '''
  A command to perform an export operation

  '''

  name = 'export'

  allowed_parents = ['integrate']

  def run(self):
    '''
    Run the export command

    '''
    from os.path import exists
    from os.path import join
    import shutil

    # Get output filename if set
    result_filename = self.phil_scope.get(diff=False).extract().mtz.hklout
    if result_filename is None:
      result_filename = "%d_integrated.mtz" % self.state.index

    # Set filenames
    self.filenames = {
      'input.experiments'  : self.state.parent.experiments,
      'input.reflections'  : self.state.parent.reflections,
      'mtz.hklout'         : join(self.state.directory, "reflections.mtz"),
      'output.log'         : join(self.state.directory, "info.log"),
      'output.debug_log'   : join(self.state.directory, "debug.log"),
    }
    self.phil_scope.set("format=mtz")
    for name, value in self.filenames.iteritems():
      self.phil_scope.set('%s=%s' % (name, value))

    # Set the parameter filename and write to file
    with open(self.state.parameters, "w") as outfile:
      outfile.write(self.phil_scope.get(diff=True).as_str())
    outfile.close()

    # Run the command
    run_external_command(['dials.export', self.state.parameters],
                         self.state.output)

    # Check the files exist
    self.check_files_exist(self.filenames.values())

    # Copy the resulting mtz file to the working directory
    shutil.copy2(self.filenames['mtz.hklout'], result_filename)


class ApplicationState(object):
  '''
  A class to hold all the application state

  '''

  # The command classes
  CommandClass = {
    'import'                             : Import,
    'find_spots'                         : FindSpots,
    'discover_better_experimental_model' : DiscoverBetterExperimentalModel,
    'index'                              : Index,
    'refine_bravais_settings'            : RefineBravaisSettings,
    'reindex'                            : Reindex,
    'refine'                             : Refine,
    'integrate'                          : Integrate,
    'export'                             : Export
  }

  class Memento(object):
    '''
    Class to init from state

    '''
    def __init__(self,
                 directory=None,
                 current=None,
                 commands=None,
                 counter=None,
                 mode=None):
      self.directory = directory
      self.current = current
      self.commands = commands
      self.counter = counter
      self.mode = mode

  def __init__(self, directory=None, memento=None):
    '''
    Initialise the state

    :param directory: The output directory

    '''

    # Get the global parameters
    self.parameters = GlobalParameterManager()

    # Init the state
    if memento is None:
      self.counter = Counter()
      self.current = InitialState()
      self.command_tree = CommandTree(self.current, self.counter)
      self.directory = directory
      self.mode = 'import'
    else:
      self.counter = Counter(memento.counter)
      self.command_tree = CommandTree(memento.commands, self.counter)
      self.current = self.command_tree.goto(memento.current)
      self.directory = memento.directory
      self.mode = memento.mode

  def run(self):
    '''
    Run the command for the given mode

    '''
    # Create the command
    command = self.CommandClass[self.mode](
      parent     = self.current,
      index      = self.counter.current(),
      phil_scope = self.parameters[self.mode],
      workspace  = self.directory)

    # Increment the counter
    self.counter.incr()

    # Apply the command
    self.current = command.apply()

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

  def as_dict(self):
    '''
    Return the application state as a dictionary

    '''

    # The application state as a dictionary
    obj = {
      'counter'   : self.counter.current(),
      'current'   : self.current.index,
      'mode'      : self.mode,
      'directory' : self.directory,
      'commands'  : self.command_tree.root.as_dict()
    }

    # Return the dictionary
    return obj

  @classmethod
  def from_dict(cls, dictionary):
    '''
    Load the class from a dictionary

    '''
    memento = cls.Memento(
      counter   = dictionary['counter'],
      current   = dictionary['current'],
      mode      = dictionary['mode'],
      directory = dictionary['directory'],
      commands  = CommandState.from_dict(dictionary['commands']))
    return cls(memento=memento)

  def dump(self, filename):
    '''
    Dump the state to file

    :param filename: The filename

    '''
    import json
    with open(filename, "w") as outfile:
      json.dump(self.as_dict(), outfile, indent=2, ensure_ascii=True)

  @classmethod
  def load(cls, filename):
    '''
    Load the state from file

    :param filename: The filename
    :return: The state object

    '''
    import json
    def _decode_list(data):
      rv = []
      for item in data:
        if isinstance(item, unicode):
          item = item.encode('utf-8')
        elif isinstance(item, list):
          item = _decode_list(item)
        elif isinstance(item, dict):
          item = _decode_dict(item)
        rv.append(item)
      return rv

    def _decode_dict(data):
      rv = {}
      for key, value in data.iteritems():
        if isinstance(key, unicode):
          key = key.encode('utf-8')
        if isinstance(value, unicode):
          value = value.encode('utf-8')
        elif isinstance(value, list):
          value = _decode_list(value)
        elif isinstance(value, dict):
          value = _decode_dict(value)
        rv[key] = value
      return rv
    with open(filename) as infile:
      return cls.from_dict(json.load(infile, object_hook=_decode_dict))


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
    if recover and exists(state_filename):
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

  def undo_parameters(self):
    '''
    Undo the last parameter changes

    '''
    self.state.parameters[self.get_mode()].undo()
    self.state.dump(self.state_filename)

  def redo_parameters(self):
    '''
    Redo the last parameter changes

    '''
    self.state.parameters[self.get_mode()].redo()
    self.state.dump(self.state_filename)

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
