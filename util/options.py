#!/usr/bin/env python
#
# options.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division
import optparse


class ConfigWriter(object):
  '''Class to write configuration to file.'''

  def __init__(self, master_phil):
    '''Initialise with the master phil.'''
    self._master_phil = master_phil

  def write(self, params, filename):
    '''Write the configuration to file.'''
    # Get the modified phil
    modified_phil = self._master_phil.format(python_object=params)

    # Get the phil as text
    text = modified_phil.as_str()

    # Write the text to file
    with open(filename, 'w') as f:
      f.write(text)


class Importer(object):
  ''' A class to import the command line arguments. '''

  def __init__(self, args,
               read_datablocks=False,
               read_experiments=False,
               read_reflections=False,
               read_datablocks_from_images=False,
               check_format=True,
               verbose=False):
    ''' Parse the arguments. Populates its instance attributes in an intelligent
    way from the arguments in args.

    If include is set, only those items set will be tried. If not, then if
    exclude is set, then those items will not be tested.

    These are the types we can import:
     - images: a list of images
     - reflections : a list of reflections
     - datablocks : a list of datablocks
     - experiments: a list of experiments

    Params:
      args The arguments to parse
      read_datablocks Try to read the datablocks
      read_experiments Try to read the experiments
      read_reflections Try to read the reflections
      read_datablocks_from_images Try to read the datablocks from images
      check_format Check the format when reading images
      verbose True/False print out some stuff

    '''

    # Initialise output
    self.datablocks = []
    self.experiments = []
    self.reflections = []
    self.unhandled = args

    # First try to read image files
    if read_datablocks_from_images:
      self.unhandled = self.try_read_datablocks_from_images(
        self.unhandled, verbose)

    # Second try to read data block files
    if read_datablocks:
      self.unhandled = self.try_read_datablocks(
        self.unhandled, check_format, verbose)

    # Third try to read experiment files
    if read_experiments:
      self.unhandled = self.try_read_experiments(
        self.unhandled, check_format, verbose)

    # Fourth try to read reflection files
    if read_reflections:
      self.unhandled = self.try_read_reflections(
        self.unhandled, verbose)

  def try_read_datablocks_from_images(self, args, verbose):
    ''' Try to import images. '''
    from dxtbx.datablock import DataBlockFactory
    from dials.phil import FilenameDataWrapper, DataBlockConverters
    unhandled = []
    datablocks = DataBlockFactory.from_filenames(
      args, verbose=verbose, unhandled=unhandled)
    if len(datablocks) > 0:
      filename = "<image files>"
      obj = FilenameDataWrapper(filename, datablocks)
      DataBlockConverters.cache[filename] = obj
      self.datablocks.append(obj)
    return unhandled

  def try_read_datablocks(self, args, check_format, verbose):
    ''' Try to import imagesets. '''
    from dials.phil import DataBlockConverters
    converter = DataBlockConverters(check_format)
    unhandled = []
    for argument in args:
      try:
        self.datablocks.append(converter.from_string(argument))
      except Exception:
        unhandled.append(argument)
    return unhandled

  def try_read_experiments(self, args, check_format, verbose):
    ''' Try to import experiments. '''
    from dials.phil import ExperimentListConverters
    converter = ExperimentListConverters(check_format)
    unhandled = []
    for argument in args:
      try:
        self.experiments.append(converter.from_string(argument))
      except Exception:
        unhandled.append(argument)
    return unhandled

  def try_read_reflections(self, args, verbose):
    ''' Try to import reflections. '''
    from dials.phil import ReflectionTableConverters
    converter = ReflectionTableConverters()
    unhandled = []
    for argument in args:
      try:
        self.reflections.append(converter.from_string(argument))
      except Exception:
        unhandled.append(argument)
    return unhandled


class PhilCommandParser(object):
  ''' A class to parse phil parameters from positional arguments '''

  def __init__(self,
               phil=None,
               read_datablocks=False,
               read_experiments=False,
               read_reflections=False,
               read_datablocks_from_images=False,
               check_format=True):
    ''' Initialise the class. '''
    from dials.phil import parse

    # Set the system phil scope
    if phil is None:
      self._system_phil = parse("")
    else:
      self._system_phil = phil

    # Set the flags
    self._read_datablocks = read_datablocks
    self._read_experiments = read_experiments
    self._read_reflections = read_reflections
    self._read_datablocks_from_images = read_datablocks_from_images
    self._check_format = check_format

    # Adopt the input scope
    input_phil_scope = self._generate_input_scope()
    if input_phil_scope is not None:
      self.system_phil.adopt_scope(input_phil_scope)

    # Set the working phil scope
    self._phil = self.system_phil.fetch(source=parse(""))

  @property
  def phil(self):
    '''Get the phil object'''
    return self._phil

  @property
  def system_phil(self):
    '''Get the system phil.'''
    return self._system_phil

  @property
  def diff_phil(self):
    ''' Get the diff phil. '''
    return self.system_phil.fetch_diff(source=self.phil)

  def parse_args(self, args, verbose=False, return_unhandled=True):
    ''' Parse the command line arguments. '''

    # Try to import everything
    importer = Importer(
      args,
      read_datablocks=self._read_datablocks,
      read_experiments=self._read_experiments,
      read_reflections=self._read_reflections,
      read_datablocks_from_images=self._read_datablocks_from_images,
      check_format=self._check_format,
      verbose=verbose)

    # Add the cached arguments
    args = importer.unhandled
    for obj in importer.datablocks:
      args.append("input.datablock=%s" % obj.filename)
    for obj in importer.experiments:
      args.append("input.experiments=%s" % obj.filename)
    for obj in importer.reflections:
      args.append("input.reflections=%s" % obj.filename)

    # Parse the command line phil parameters
    interpretor = self.system_phil.command_line_argument_interpreter()
    if return_unhandled is True:
      self._phil, args = interpretor.process_and_fetch(args,
        custom_processor="collect_remaining")
    else:
      processed = interpretor.process(args=args)
      self._phil, unused = self._system_phil.fetch(
        sources=processed,
        track_unused_definitions=True)
      if (len(unused) != 0):
        raise RuntimeError((
          '\n'
          ' The following phil parameters were not recognised:\n'
          '  %s\n'
        ) % '\n  '.join(map(str,unused)))
      args = []
    return self.phil.extract(), args

  def _generate_input_scope(self):
    ''' Generate the required input scope. '''
    from dials.phil import parse

    # Create the input scope
    require_input_scope = (
      self._read_datablocks or
      self._read_experiments or
      self._read_reflections or
      self._read_datablocks_from_images)
    if not require_input_scope:
      return None
    input_phil_scope = parse('input {}')
    main_scope = input_phil_scope.get_without_substitution("input")
    assert(len(main_scope) == 1)
    main_scope = main_scope[0]

    # Add the datablock phil scope
    if self._read_datablocks or self._read_datablocks_from_images:
      phil_scope = parse('''
        datablock = None
          .type = datablock(check_format=%r)
          .multiple = True
          .help = "The datablock file path"
      ''' % self._check_format)
      main_scope.adopt_scope(phil_scope)

    # Add the experiments phil scope
    if self._read_experiments:
      phil_scope = parse('''
        experiments = None
          .type = experiment_list(check_format=%r)
          .multiple = True
          .help = "The experiment list file path"
      ''' % self._check_format)
      main_scope.adopt_scope(phil_scope)

    # Add the reflections scope
    if self._read_reflections:
      phil_scope = parse('''
        reflections = None
          .type = reflection_table
          .multiple = True
          .help = "The reflection table file path"
      ''')
      main_scope.adopt_scope(phil_scope)

    # Return the input scope
    return input_phil_scope


class OptionParserBase(optparse.OptionParser, object):
  ''' The base class for the option parser. '''

  def __init__(self,
               config_options=False,
               sort_options=False,
               **kwargs):
    '''Initialise the class.'''

    # Initialise the option parser
    super(OptionParserBase, self).__init__(**kwargs)

    # Add an option to show configuration parameters
    if config_options:
      self.add_option(
        '-c', '--show-config',
        action='store_true',
        default=False,
        dest='show_config',
        help='Show the configuration parameters.')
      self.add_option(
        '-a', '--attributes-level',
        default=2,
        type='int',
        dest='attributes_level',
        help='Set the attributes level for showing configuration parameters')
      self.add_option(
        '-e', '--expert-level',
        type='int',
        default=0,
        dest='expert_level',
        help='Set the expert level for showing configuration parameters')

    # Add an option to sort
    if sort_options:
      self.add_option(
        '-s', '--sort',
        action='store_true',
        dest='sort',
        default=False,
        help='Sort the arguments')

    # Set a verbosity parameter
    self.add_option(
      '-v',
      action='count',
      default=0,
      dest='verbose',
      help='Set the verbosity')

  def parse_args(self, args=None):
    '''Parse the command line arguments and get system configuration.'''
    import sys
    import select

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = super(OptionParserBase, self).parse_args(args=args)

    # Read stdin if data is available
    r, w, x = select.select([sys.stdin], [], [], 0)
    if len(r) > 0:
      args.extend([l.strip() for rr in r for l in rr.readlines()])

    # Maybe sort the data
    if hasattr(options, "sort") and options.sort:
      args = sorted(args)

    # Return the parameters
    return options, args

  def format_epilog(self, formatter):
    ''' Don't do formatting on epilog. '''
    if self.epilog is None:
      return ''
    return self.epilog


class OptionParser(OptionParserBase):
  '''A class to parse command line options and get the system configuration.
  The class extends optparse.OptionParser to include the reading of phil
  parameters.'''

  def __init__(self,
               phil=None,
               read_datablocks=False,
               read_experiments=False,
               read_reflections=False,
               read_datablocks_from_images=False,
               check_format=True,
               sort_options=False,
               **kwargs):
    '''Initialise the class.'''
    from dials.phil import parse

    # Create the phil parser
    self._phil_parser = PhilCommandParser(
      phil=phil,
      read_datablocks=read_datablocks,
      read_experiments=read_experiments,
      read_reflections=read_reflections,
      read_datablocks_from_images=read_datablocks_from_images,
      check_format=check_format)

    # Initialise the option parser
    super(OptionParser, self).__init__(
      sort_options=sort_options,
      config_options=self.system_phil.as_str() != '',
      **kwargs)

  def parse_args(self, args=None, show_diff_phil=False, return_unhandled=False):
    '''Parse the command line arguments and get system configuration.'''
    import sys

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = super(OptionParser, self).parse_args(args=args)

    # Show config
    if hasattr(options, 'show_config') and options.show_config:
      print (
        'Showing configuration parameters with:\n'
        '  attributes_level = %d\n'
        '  expert_level = %d\n' % (
          options.attributes_level,
          options.expert_level))
      print self.phil.as_str(
        expert_level=options.expert_level,
        attributes_level=options.attributes_level)
      exit(0)

    # Parse the phil parameters
    params, args = self._phil_parser.parse_args(
      args, options.verbose > 0,
      return_unhandled=return_unhandled)

    # Print the diff phil
    if show_diff_phil:
      diff_phil_str = self.diff_phil.as_str()
      if (diff_phil_str is not ''):
        print 'The following parameters have been modified:\n'
        print diff_phil_str

    # Return the parameters
    if return_unhandled:
      return params, options, args
    else:
      if args:
        raise RuntimeError((
          'This should not be reached!\n'
          'Unhandled arguments: %s') % (' '.join(args)))
    return params, options

  @property
  def phil(self):
    '''Get the phil object'''
    return self._phil_parser.phil

  @property
  def system_phil(self):
    '''Get the system phil.'''
    return self._phil_parser.system_phil

  @property
  def diff_phil(self):
    ''' Get the diff phil. '''
    return self._phil_parser.diff_phil

  def _strip_rst_markup(self, text):
    return text.replace("::", ":")

  def format_help(self, formatter=None):
    result = super(OptionParser, self).format_help(formatter=formatter)
    return self._strip_rst_markup(result)


def flatten_reflections(filename_object_list):
  result = []
  for i in range(len(filename_object_list)):
    result.append(filename_object_list[i].data)
  return result

def flatten_datablocks(filename_object_list):
  result = []
  for i in range(len(filename_object_list)):
    result.extend(filename_object_list[i].data)
  return result

def flatten_experiments(filename_object_list):
  from dxtbx.model.experiment.experiment_list import ExperimentList
  result = ExperimentList()
  for i in range(len(filename_object_list)):
    result.extend(filename_object_list[i].data)
  return result
