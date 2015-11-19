from __future__ import division
from cmd import Cmd


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
      input {
        datablock = None
          .type = str
      }
      include scope dials.command_line.find_spots.phil_scope
    ''', process_includes=True)
    super(FindSpotsParameterManager, self).__init__(phil_scope)


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
      input {
        datablock = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.index.phil_scope
    ''', process_includes=True)
    super(IndexParameterManager, self).__init__(phil_scope)


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
      input {
        experiments = None
          .type = str
        reflections = None
          .type = str
      }
      include scope dials.command_line.export_mtz.phil_scope
    ''', process_includes=True)
    super(ExportParameterManager, self).__init__(phil_scope)


class GlobalParameterManager(object):
  '''
  Class to hold all parameter managers

  '''

  def __init__(self):
    '''
    Init everything

    '''
    self.import_parameters = ImportParameterManager()
    self.find_spots_parameters = FindSpotsParameterManager()
    self.index_parameters = IndexParameterManager()
    self.refine_parameters = RefineParameterManager()
    self.integrate_parameters = IntegrateParameterManager()
    self.export_parameters = ExportParameterManager()



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

    # Init the result
    self.result = None

    # Save the info
    self.action = action
    self.parameters = copy.deepcopy(parameters)
    if directory is not None:
      self.directory = join(directory, "%d_%s" % (self.index, self.action))
    else:
      self.directory = None

  def __iter__(self):
    '''
    Iterate through the children and their children

    '''
    yield self, 0
    for child in self.children:
      for node, depth in child:
        yield node, depth+1

  def done(self):
    '''
    Check if we've done the command

    :return: True/False are we done
    '''
    return self.result != None

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

    # Initialise running the command
    self.initialize()

    # Make the directory to store output
    makedirs(self.directory)

    # Set the parameter filename and write to file
    parameters = join(self.directory, "parameters.phil")
    with open(parameters, "w") as outfile:
      outfile.write(self.parameters.get(diff=True).as_str())
    outfile.close()

    # Set the output filename
    output = join(self.directory, "output.txt")

    # Run the command (override this method)
    self.run(parameters, output)

    # Grab the result (override this method)
    self.finalize()

    # Return success/failure
    return True#self.result.success


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
        buf.write(prefix[:-3])
        buf.write('  +--')
      buf.write(node.action)
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

  def apply(self):
    raise RuntimeError("Nothing to do")


class ImportImagesCommand(CommandNode):

  parent_actions = ['clean']

  def __init__(self, parent, parameters, directory):

    super(ImportImagesCommand, self).__init__(
      parent, 'import', parameters, directory)

  def initialize(self):
    from os.path import join
    self.output_datablock_filename = join(self.directory, "datablock.json")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.parameters.set("output=%s" % self.output_datablock_filename)
    self.parameters.set("log=%s" % self.output_info_log_filename)
    self.parameters.set("debug_log=%s" % self.output_debug_log_filename)

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running import: for output see %s" % output
    easy_run.fully_buffered([
      'dials.import',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    for filename in [
        self.output_datablock_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)


class FindSpotsCommand(CommandNode):

  parent_actions = ['import']

  def __init__(self, parent, parameters, directory):
    super(FindSpotsCommand, self).__init__(
      parent, 'find_spots', parameters, directory)

  def initialize(self):
    from os.path import join
    self.input_datablock_filename = self.parent.output_datablock_filename
    self.output_reflections_filename = join(self.directory, "reflections.pickle")
    self.output_datablock_filename = join(self.directory, "datablock.json")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.parameters.set("input.datablock=%s" % self.input_datablock_filename)
    self.parameters.set("output.reflections=%s" % self.output_reflections_filename)
    self.parameters.set("output.datablock=%s" % self.output_datablock_filename)
    self.parameters.set("output.log=%s" % self.output_info_log_filename)
    self.parameters.set("output.debug_log=%s" % self.output_debug_log_filename)

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running find_spots: for output see %s" % output
    easy_run.fully_buffered([
      'dials.find_spots',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    for filename in [
        self.output_reflections_filename,
        self.output_datablock_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)


class IndexCommand(CommandNode):

  parent_actions = ['find_spots']

  def __init__(self, parent, parameters, directory):
    super(IndexCommand, self).__init__(
      parent, 'index', parameters, directory)

  def initialize(self):
    from os.path import join
    self.input_datablock_filename = self.parent.output_datablock_filename
    self.input_reflections_filename = self.parent.output_reflections_filename
    self.output_reflections_filename = join(self.directory, "reflections.pickle")
    self.output_experiments_filename = join(self.directory, "experiments.json")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.parameters.set("input.datablock=%s" % self.input_datablock_filename)
    self.parameters.set("input.reflections=%s" % self.input_reflections_filename)
    self.parameters.set("output.reflections=%s" % self.output_reflections_filename)
    self.parameters.set("output.experiments=%s" % self.output_experiments_filename)
    self.parameters.set("output.log=%s" % self.output_info_log_filename)
    self.parameters.set("output.debug_log=%s" % self.output_debug_log_filename)

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running index: for output see %s" % output
    easy_run.fully_buffered([
      'dials.index',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    for filename in [
        self.output_reflections_filename,
        self.output_experiments_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)


class RefineCommand(CommandNode):

  parent_actions = ['index', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    super(RefineCommand, self).__init__(
      parent, 'refine', parameters, directory)

  def initialize(self):
    from os.path import join
    self.input_experiments_filename = self.parent.output_experiments_filename
    self.input_reflections_filename = self.parent.output_reflections_filename
    self.output_reflections_filename = join(self.directory, "reflections.pickle")
    self.output_experiments_filename = join(self.directory, "experiments.json")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.output_matches_filename = join(self.directory, "matches.pickle")
    self.output_centroids_filename = join(self.directory, "centroids.txt")
    # self.output_parameter_table_filename = join(self.directory, "parameter_table.txt")
    self.output_history_filename = join(self.directory, "history.txt")
    self.parameters.set("input.experiments=%s" % self.input_experiments_filename)
    self.parameters.set("input.reflections=%s" % self.input_reflections_filename)
    self.parameters.set("output.reflections=%s" % self.output_reflections_filename)
    self.parameters.set("output.experiments=%s" % self.output_experiments_filename)
    self.parameters.set("output.log=%s" % self.output_info_log_filename)
    self.parameters.set("output.debug_log=%s" % self.output_debug_log_filename)
    self.parameters.set("output.matches=%s" % self.output_matches_filename)
    self.parameters.set("output.centroids=%s" % self.output_centroids_filename)
    # self.parameters.set("output.parameter_table=%s" % self.output_parameter_table_filename)
    self.parameters.set("output.history=%s" % self.output_history_filename)

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running refine: for output see %s" % output
    easy_run.fully_buffered([
      'dials.refine',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    for filename in [
        self.output_reflections_filename,
        self.output_experiments_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename,
        self.output_matches_filename,
        self.output_centroids_filename,
        # self.output_parameter_table_filename,
        self.output_history_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)


class IntegrateCommand(CommandNode):

  parent_actions = ['index', 'refine', 'integrate']

  def __init__(self, parent, parameters, directory):
    super(IntegrateCommand, self).__init__(
      parent, 'integrate', parameters, directory)

  def initialize(self):
    from os.path import join
    self.input_experiments_filename = self.parent.output_experiments_filename
    self.input_reflections_filename = self.parent.output_reflections_filename
    self.output_reflections_filename = join(self.directory, "reflections.pickle")
    self.output_experiments_filename = join(self.directory, "experiments.json")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.output_report_filename = join(self.directory, "report.json")
    self.parameters.set("input.experiments=%s" % self.input_experiments_filename)
    self.parameters.set("input.reflections=%s" % self.input_reflections_filename)
    self.parameters.set("output.reflections=%s" % self.output_reflections_filename)
    self.parameters.set("output.experiments=%s" % self.output_experiments_filename)
    self.parameters.set("output.log=%s" % self.output_info_log_filename)
    self.parameters.set("output.debug_log=%s" % self.output_debug_log_filename)
    self.parameters.set("output.report=%s" % self.output_report_filename)
    self.parameters.set("output.phil=None")

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running integrate: for output see %s" % output
    easy_run.fully_buffered([
      'dials.integrate',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    for filename in [
        self.output_reflections_filename,
        self.output_experiments_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename,
        self.output_report_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)


class ExportCommand(CommandNode):

  parent_actions = ['integrate']

  def __init__(self, parent, parameters, directory):
    super(ExportCommand, self).__init__(
      parent, 'export', parameters, directory)

  def initialize(self):
    from os.path import join
    self.input_experiments_filename = self.parent.output_experiments_filename
    self.input_reflections_filename = self.parent.output_reflections_filename
    self.output_reflections_filename = join(self.directory, "reflections.mtz")
    self.output_info_log_filename = join(self.directory, "info.log")
    self.output_debug_log_filename = join(self.directory, "debug.log")
    self.parameters.set("input.experiments=%s" % self.input_experiments_filename)
    self.parameters.set("input.reflections=%s" % self.input_reflections_filename)
    self.parameters.set("hklout=%s" % self.output_reflections_filename)
    self.parameters.set("log=%s" % self.output_info_log_filename)
    self.parameters.set("debug_log=%s" % self.output_debug_log_filename)

  def run(self, parameters, output):
    from libtbx import easy_run
    print "Running export: for output see %s" % output
    easy_run.fully_buffered([
      'dials.export_mtz',
      parameters,
      '>',
      output])#.raise_if_errors()

  def finalize(self):
    from os.path import exists
    import shutil
    for filename in [
        self.output_reflections_filename,
        self.output_info_log_filename,
        self.output_debug_log_filename]:
      if not exists(filename):
        raise RuntimeError("File %s could not be found" % filename)

    # Copy the resulting mtz file to the working directory
    result_filename = "%d_integrated.mtz" % self.index
    shutil.copy2(self.output_reflections_filename, result_filename)


class CommandRunner(object):

  def __init__(self, directory):
    # Create the parameters
    self.parameters = GlobalParameterManager()

    # Set the initial state to current
    self.current = InitialState()

    # Create the command tree
    self.command_tree = CommandTree(self.current)

    # Save the parameters and directory
    self.directory = directory

  def run_command(self, parameters, Class):
    self.current = Class(
      self.current,
      parameters,
      self.directory)
    self.current.apply()
    return self.current.index

  def import_images(self):
    return self.run_command(
      self.parameters.import_parameters,
      ImportImagesCommand)

  def find_spots(self):
    return self.run_command(
      self.parameters.find_spots_parameters,
      FindSpotsCommand)

  def index(self):
    return self.run_command(
      self.parameters.index_parameters,
      IndexCommand)

  def refine(self):
    return self.run_command(
      self.parameters.refine_parameters,
      RefineCommand)

  def integrate(self):
    return self.run_command(
      self.parameters.integrate_parameters,
      IntegrateCommand)

  def export(self):
    return self.run_command(
      self.parameters.export_parameters,
      ExportCommand)

  def goto(self, index):
    self.current = self.command_tree.goto(index)

  def string(self):
    return self.command_tree.string(current=self.current.index)

  def dump(self, filename):
    import cPickle as pickle
    with open(filename, "w") as outfile:
      pickle.dump(self, outfile)

  @classmethod
  def load(Class, filename):
    import cPickle as pickle
    with open(filename) as infile:
      return pickle.load(infile)


class Controller(object):

  mode_list = ['import', 'find_spots', 'index', 'refine', 'integrate', 'export']

  def __init__(self, directory="output", state_filename="dials.state", recover=True):
    from os.path import exists, abspath

    # Set some stuff
    self.state_filename = state_filename
    self.mode = "import"

    # Read state if available
    if recover == True and exists(state_filename):
      self.state = CommandRunner.load(state_filename)
      print "Recovered state from %s" % state_filename
      self.show()
    else:
      self.state = CommandRunner(abspath(directory))

  def program(self, program):
    if program not in self.mode_list:
      raise RuntimeError('Unknown mode: %s' % program)
    self.mode = program

  def set(self, parameters, short_syntax=False):
    from libtbx.utils import Sorry
    try:
      self.parameters().set(parameters, short_syntax=short_syntax)
      self.state.dump(self.state_filename)
    except Sorry, e:
      print e
      return False
    return True

  def reset(self):
    self.parameters().reset()
    self.state.dump(self.state_filename)
    return True

  def get(self, diff=True):
    return self.parameters().get(diff=diff)

  def show(self):
    print self.state.string()

  def goto(self, index):
    self.state.goto(index)
    self.show()

  def parameters(self):
    if self.mode == 'import':
      return self.state.parameters.import_parameters
    elif self.mode == 'find_spots':
      return self.state.parameters.find_spots_parameters
    elif self.mode == 'index':
      return self.state.parameters.index_parameters
    elif self.mode == 'refine':
      return self.state.parameters.refine_parameters
    elif self.mode == 'integrate':
      return self.state.parameters.integrate_parameters
    elif self.mode == 'export':
      return self.state.parameters.export_parameters
    else:
      raise RuntimeError('Unknown mode: %s' % self.mode)

  def run(self):
    if self.mode == 'import':
      result = self.state.import_images()
    elif self.mode == 'find_spots':
      result = self.state.find_spots()
    elif self.mode == 'index':
      result = self.state.index()
    elif self.mode == 'refine':
      result = self.state.refine()
    elif self.mode == 'integrate':
      result = self.state.integrate()
    elif self.mode == 'export':
      result = self.state.export()
    else:
      raise RuntimeError('Unknown mode: %s' % self.mode)
    self.state.dump(self.state_filename)
    return result


class Console(Cmd):
  ''' Interactive console. '''

  intro = 'DIALS interactive mode (type "help" for documentation)'
  prompt = ">> "

  def __init__(self):

    Cmd.__init__(self)

    self.controller = Controller()
    self.prompt = "%s >> " % self.controller.mode

  def emptyline(self):
    ''' Do nothing on empty line '''
    pass

  def default(self, line):
    try:
      self.do_set(line)
    except Exception:
      return Cmd.default(line)

  def do_mode(self, mode):
    ''' Set the program mode '''
    try:
      self.controller.program(mode)
      self.prompt = "%s >> " % self.controller.mode
    except Exception, e:
      print e

  def complete_mode(self, text, line, begidx, endidx):
    return [i for i in self.controller.mode_list if i.startswith(text)]

  def do_reset(self, line):
    ''' Reset parameters to default. '''
    self.controller.reset()

  def do_set(self, parameter):
    ''' Set a phil parameter '''
    self.controller.set(parameter, short_syntax=True)

  def do_load(self, filename):
    ''' Load a phil parameter file '''
    with open(filename) as infile:
      self.controller.set(infile.read())

  def do_get(self, line):
    ''' Show the modified parameters. '''
    print self.controller.get(diff=True).as_str()

  def do_all(self, expert_level):
    ''' Show all the possible parameters '''
    print self.controller.get(diff=False).as_str(expert_level=str(expert_level))

  def do_run(self, line):
    ''' Run a program '''
    self.controller.run()

  def do_goto(self, line):
    ''' Goto a particular history state '''
    self.controller.goto(int(line))

  def do_history(self, line):
    ''' Show the history. '''
    self.controller.show()

  def do_import(self, params):
    ''' Imperative import command '''
    self.run_as_imperative("import", params)

  def do_find_spots(self, params):
    ''' Imperative find_spots command '''
    self.run_as_imperative("find_spots", params)

  def do_index(self, params):
    ''' Imperative index command '''
    self.run_as_imperative("index", params)

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

  def run_as_imperative(self, mode, params):
    self.do_mode(mode)
    self.do_set(params)
    self.do_run("")


if __name__ == '__main__':
  console = Console()
  console.cmdloop()
