
#
# options.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division
import optparse

from libtbx.utils import Sorry

import libtbx.phil
tolerance_phil_scope = libtbx.phil.parse('''
tolerance
    .help = "Tolerances used to determine shared models"
    .expert_level = 2
  {

  beam {

    wavelength = 1e-6
      .type = float(value_min=0.0)
      .help = "The wavelength tolerance"

    direction = 1e-6
      .type = float(value_min=0.0)
      .help = "The direction tolerance"

    polarization_normal = 1e-6
      .type = float(value_min=0.0)
      .help = "The polarization normal tolerance"

    polarization_fraction = 1e-6
      .type = float(value_min=0.0)
      .help = "The polarization fraction tolerance"

  }

  detector {

    fast_axis = 1e-6
      .type = float(value_min=0.0)
      .help = "The fast axis tolerance"

    slow_axis = 1e-6
      .type = float(value_min=0.0)
      .help = "The slow axis tolerance"

    origin = 5e-2
      .type = float(value_min=0.0)
      .help = "The origin tolerance"

  }

  goniometer {

    rotation_axis = 1e-6
      .type = float(value_min=0.0)
      .help = "The rotation axis tolerance"

    fixed_rotation = 1e-6
      .type = float(value_min=0.0)
      .help = "The fixed rotation tolerance"

    setting_rotation = 1e-6
      .type = float(value_min=0.0)
      .help = "The setting rotation tolerance"

  }

  scan {

    oscillation = 0.01
      .type = float(value_min=0.0)
      .help = "The oscillation tolerance for the scan"

  }
}
''')

geometry_phil_scope = libtbx.phil.parse('''
geometry
  .help = "Allow overrides of experimental geometry"
  .expert_level = 2
{
  include scope dxtbx.model.beam.beam_phil_scope
  include scope dxtbx.model.detector.detector_phil_scope
  include scope dxtbx.model.goniometer.goniometer_phil_scope
  include scope dxtbx.model.scan.scan_phil_scope

  convert_stills_to_sweeps = False
    .type = bool
    .help = "When overriding the scan, convert stills into sweeps"
    .short_caption = "Convert stills into sweeps"

  convert_sweeps_to_stills = False
    .type = bool
    .help = "When overriding the scan, convert sweeps into stills"
    .short_caption = "Convert sweeps into stills"
}
''', process_includes=True)


format_phil_scope = libtbx.phil.parse('''
format
  .help = "Options to pass to the Format class"
  .expert_level = 2
{
  dynamic_shadowing = False
    .type = bool
    .help = "Enable dynamic shadowing"
  multi_panel = False
    .type = bool
    .help = "Enable a multi-panel detector model."
            "(Not supported by all detector formats)"
}
''')


class ConfigWriter(object):
  '''Class to write configuration to file.'''

  def __init__(self, master_phil):
    '''
    Initialise with the master phil.

    :param master_phil: The master phil scope

    '''
    self._master_phil = master_phil

  def write(self, params, filename):
    '''
    Write the configuration to file.

    :param params: The input phil parameters
    :param filename: The output filename

    '''
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
               verbose=False,
               compare_beam=None,
               compare_detector=None,
               compare_goniometer=None,
               scan_tolerance=None,
               format_kwargs=None):
    '''
    Parse the arguments. Populates its instance attributes in an intelligent way
    from the arguments in args.

    If include is set, only those items set will be tried. If not, then if
    exclude is set, then those items will not be tested.

    These are the types we can import:
     - images: a list of images
     - reflections : a list of reflections
     - datablocks : a list of datablocks
     - experiments: a list of experiments

    :param args: The arguments to parse
    :param read_datablocks: Try to read the datablocks
    :param read_experiments: Try to read the experiments
    :param read_reflections: Try to read the reflections
    :param read_datablocks_from_images: Try to read the datablocks from images
    :param check_format: Check the format when reading images
    :param verbose: True/False print out some stuff

    '''

    # Initialise output
    self.datablocks = []
    self.experiments = []
    self.reflections = []
    self.unhandled = args

    # First try to read image files
    if read_datablocks_from_images:
      self.unhandled = self.try_read_datablocks_from_images(
        self.unhandled,
        verbose,
        compare_beam,
        compare_detector,
        compare_goniometer,
        scan_tolerance,
        format_kwargs)

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

  def try_read_datablocks_from_images(self,
                                      args,
                                      verbose,
                                      compare_beam,
                                      compare_detector,
                                      compare_goniometer,
                                      scan_tolerance,
                                      format_kwargs):
    '''
    Try to import images.

    :param args: The input arguments
    :param verbose: Print verbose output
    :return: Unhandled arguments

    '''
    from dxtbx.datablock import DataBlockFactory
    from dials.phil import FilenameDataWrapper, DataBlockConverters
    from glob import glob

    # If filenames contain wildcards, expand
    args_new = []
    for arg in args:
      if "*" in arg:
        args_new.extend(glob(arg))
      else:
        args_new.append(arg)
    args = args_new

    unhandled = []
    datablocks = DataBlockFactory.from_filenames(
      args,
      verbose=verbose,
      unhandled=unhandled,
      compare_beam=compare_beam,
      compare_detector=compare_detector,
      compare_goniometer=compare_goniometer,
      scan_tolerance=scan_tolerance,
      format_kwargs=format_kwargs)
    if len(datablocks) > 0:
      filename = "<image files>"
      obj = FilenameDataWrapper(filename, datablocks)
      DataBlockConverters.cache[filename] = obj
      self.datablocks.append(obj)
    return unhandled

  def try_read_datablocks(self, args, check_format, verbose):
    '''
    Try to import imagesets.

    :param args: The input arguments
    :param check_format: True/False check the image format
    :param verbose: Print verbose output
    :returns: Unhandled arguments

    '''
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
    '''
    Try to import experiments.

    :param args: The input arguments
    :param check_format: True/False check the image format
    :param verbose: Print verbose output
    :returns: Unhandled arguments

    '''
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
    ''' Try to import reflections.

    :param args: The input arguments
    :param verbose: Print verbose output
    :returns: Unhandled arguments

    '''
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
    '''
    Initialise the parser.

    :param phil: The phil scope
    :param read_datablocks: Try to read the datablocks
    :param read_experiments: Try to read the experiments
    :param read_reflections: Try to read the reflections
    :param read_datablocks_from_images: Try to read the datablocks from images
    :param check_format: Check the format when reading images

    '''
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
    '''
    Get the phil object

    :return: The phil scope

    '''
    return self._phil

  @property
  def system_phil(self):
    '''
    Get the system phil.

    :return: The system phil scope

    '''
    return self._system_phil

  @property
  def diff_phil(self):
    '''
    Get the diff phil.

    :return: The difference phil scope

    '''
    return self.system_phil.fetch_diff(source=self.phil)

  def parse_args(self, args, verbose=False, return_unhandled=False, quick_parse=False):
    '''
    Parse the command line arguments.

    :param args: The input arguments
    :param verbose: Print verbose output
    :param return_unhandled: True/False also return unhandled arguments
    :param quick_parse: Return as fast as possible and without reading any data,
                        ignoring class constructor options.
    :return: The options and parameters and (optionally) unhandled arguments

    '''
    import os.path
    from dxtbx.datablock import BeamComparison
    from dxtbx.datablock import DetectorComparison
    from dxtbx.datablock import GoniometerComparison
    from dials.phil import parse

    # Parse the command line phil parameters
    user_phils = []
    unhandled = []
    interpretor = self.system_phil.command_line_argument_interpreter()
    for arg in args:
      if os.path.isfile(arg) and os.path.getsize(arg) > 0:
        name, ext = os.path.splitext(arg)
        if ext in ['.phil', '.param', '.params', '.eff', '.def']:
          try:
            user_phils.append(parse(file_name=arg))
          except Exception:
            if return_unhandled:
              unhandled.append(arg)
            else:
              raise
        else:
          unhandled.append(arg)
      elif arg.find("=") >= 0:
        try:
          user_phils.append(interpretor.process_arg(arg=arg))
        except Exception:
          if return_unhandled:
            unhandled.append(arg)
          else:
            raise
      else:
        unhandled.append(arg)

    # Fetch the phil parameters
    self._phil, unused = self.system_phil.fetch(
      sources=user_phils,
      track_unused_definitions=True)

    # Print if bad definitions
    if len(unused) > 0:
      msg = [item.object.as_str().strip() for item in unused]
      msg = '\n'.join(['  %s' % line for line in msg])
      raise RuntimeError('The following definitions were not recognised\n%s' % msg)

    # Extract the parameters
    params = self._phil.extract()

    # Stop at this point if quick_parse is set. A second pass may be needed.
    if quick_parse:
      return params, unhandled

    # Create some comparison functions
    if self._read_datablocks_from_images:
      compare_beam = BeamComparison(
        wavelength_tolerance=params.input.tolerance.beam.wavelength,
        direction_tolerance=params.input.tolerance.beam.direction,
        polarization_normal_tolerance=params.input.tolerance.beam.polarization_normal,
        polarization_fraction_tolerance=params.input.tolerance.beam.polarization_fraction)
      compare_detector = DetectorComparison(
        fast_axis_tolerance=params.input.tolerance.detector.fast_axis,
        slow_axis_tolerance=params.input.tolerance.detector.slow_axis,
        origin_tolerance=params.input.tolerance.detector.origin)
      compare_goniometer = GoniometerComparison(
        rotation_axis_tolerance=params.input.tolerance.goniometer.rotation_axis,
        fixed_rotation_tolerance=params.input.tolerance.goniometer.fixed_rotation,
        setting_rotation_tolerance=params.input.tolerance.goniometer.setting_rotation)
      scan_tolerance = params.input.tolerance.scan.oscillation

      # FIXME Should probably make this smarter since it requires editing here
      # and in dials.import phil scope
      try:
        format_kwargs = {
          'dynamic_shadowing' : params.format.dynamic_shadowing,
          'multi_panel' : params.format.multi_panel,
        }
      except Exception:
        format_kwargs = None
    else:
      compare_beam = None
      compare_detector = None
      compare_goniometer = None
      scan_tolerance = None
      format_kwargs = None

    # Try to import everything
    importer = Importer(
      unhandled,
      read_datablocks=self._read_datablocks,
      read_experiments=self._read_experiments,
      read_reflections=self._read_reflections,
      read_datablocks_from_images=self._read_datablocks_from_images,
      check_format=self._check_format,
      verbose=verbose,
      compare_beam=compare_beam,
      compare_detector=compare_detector,
      compare_goniometer=compare_goniometer,
      scan_tolerance=scan_tolerance,
      format_kwargs=format_kwargs)

    # Add the cached arguments
    for obj in importer.datablocks:
      params.input.datablock.append(obj)
    for obj in importer.experiments:
      params.input.experiments.append(obj)
    for obj in importer.reflections:
      params.input.reflections.append(obj)

    # Convert to phil
    self._phil = self.system_phil.format(python_object=params)

    return params, importer.unhandled

  def _generate_input_scope(self):
    '''
    Generate the required input scope.

    :return: The input phil scope

    '''
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

    # If reading images, add some more parameters
    if self._read_datablocks_from_images:
      main_scope.adopt_scope(tolerance_phil_scope)

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
    '''
    Initialise the class.

    :param config_options: True/False show configuration options
    :param sort_options: True/False show argument sorting options

    '''

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
        default=0,
        type='int',
        dest='attributes_level',
        help='Set the attributes level for showing configuration parameters')
      self.add_option(
        '-e', '--expert-level',
        type='int',
        default=0,
        dest='expert_level',
        help='Set the expert level for showing configuration parameters')
      self.add_option(
        '--export-autocomplete-hints',
        action='store_true',
        default=False,
        dest='export_autocomplete_hints',
        help=optparse.SUPPRESS_HELP)

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
      help='Increase verbosity')

  def parse_args(self, args=None, quick_parse=False):
    '''
    Parse the command line arguments and get system configuration.

    :param args: The arguments to parse.
    :returns: The options and phil parameters

    '''
    import sys
    import select
    import os

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = super(OptionParserBase, self).parse_args(args=args)

    # Read stdin if data is available
    try:
      if not quick_parse and not sys.stdin.isatty():
        args.extend(l.strip() for l in sys.stdin.readlines())
    except IOError as e:
      if e.errno != 9:
        raise # Ignore 'bad file descriptor' errors, which may be caused by nohup et al.

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
    '''
    Initialise the class.

    :param phil: The phil scope
    :param read_datablocks: Try to read the datablocks
    :param read_experiments: Try to read the experiments
    :param read_reflections: Try to read the reflections
    :param read_datablocks_from_images: Try to read the datablocks from images
    :param check_format: Check the format when reading images
    :param sort_options: Show argument sorting options

    '''
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

  def parse_args(self,
                 args=None,
                 show_diff_phil=False,
                 return_unhandled=False,
                 ignore_unhandled=False,
                 quick_parse=False):
    '''
    Parse the command line arguments and get system configuration.

    :param args: The input arguments
    :param show_diff_phil: True/False Print the diff phil
    :param return_unhandled: True/False return unhandled arguments
    :param ignore_unhandled: True/False ignore unhandled arguments
                              if return_unhandled is False
    :param quick_parse: Return as fast as possible and without reading any data,
                        ignoring class constructor options
    :return: The options and phil parameters

    '''
    import sys

    # Parse the command line arguments, this will separate out
    # options (e.g. -o, --option) and positional arguments, in
    # which phil options will be included.
    options, args = super(OptionParser, self).parse_args(
      args=args, quick_parse=quick_parse)

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

    if hasattr(options, 'export_autocomplete_hints') and options.export_autocomplete_hints:
      self._export_autocomplete_hints()
      exit(0)

    # Parse the phil parameters
    params, args = self._phil_parser.parse_args(
      args,
      options.verbose > 0,
      return_unhandled=return_unhandled,
      quick_parse=quick_parse)

    # Print the diff phil
    if show_diff_phil:
      diff_phil_str = self.diff_phil.as_str()
      if diff_phil_str is not '':
        print 'The following parameters have been modified:\n'
        print diff_phil_str

    # Return the parameters
    if return_unhandled:
      return params, options, args
    elif len(args) > 0 and not quick_parse:
      msg = 'Unable to handle the following arguments:\n'
      msg += '\n'.join(['  %s' % a for a in args])
      if any(a.endswith(".json") for a in args):
        msg += "\n"
        msg += "\nPlease check that your image data is accessible"
      if ignore_unhandled:
        print msg
      else:
        raise Sorry(msg)
    return params, options

  @property
  def phil(self):
    '''
    Get the phil object

    :returns: The phil scope

    '''
    return self._phil_parser.phil

  @property
  def system_phil(self):
    '''
    Get the system phil.

    :returns: The system phil scope

    '''
    return self._phil_parser.system_phil

  @property
  def diff_phil(self):
    '''
    Get the diff phil.

    :returns: The diff phil scope

    '''
    return self._phil_parser.diff_phil

  def _strip_rst_markup(self, text):
    '''
    Strip rst markup

    :param text: The text to strip
    :return: The stripped text

    '''
    return text.replace("::", ":")

  def format_help(self, formatter=None):
    '''
    Format the help string

    :param formatter: The formatter to use
    :return: The formatted help text

    '''
    result = super(OptionParser, self).format_help(formatter=formatter)
    return self._strip_rst_markup(result)

  def _export_autocomplete_hints(self):
    # complete list of all parameters
    parameter_list = []
    # short name -> full name expansion for unique names
    parameter_expansion_list = {}
    # full name -> list of flags for choice parameters
    parameter_choice_list = {}

    for d in self.phil.all_definitions():
      # Create complete list of all parameters
      parameter_list.append(d.path)

      # Create expansion (alias) list for unique names that are not expert commands
      if d.object.name not in parameter_expansion_list and \
         d.object.expert_level is None and \
         d.parent.expert_level is None:
        parameter_expansion_list[d.object.name] = d.path
      else:
        parameter_expansion_list[d.object.name] = None

      # Extract parameter choice lists
      if d.object.type.phil_type == "choice":
        parameter_choice_list[d.path] = \
          [w[1:] if w.startswith("*") else w for w in (str(x) for x in d.object.words) ]
      elif d.object.type.phil_type == "bool":
        parameter_choice_list[d.path] = ["true", "false"]

    def construct_completion_tree(paths):
      """ Construct a tree of parameters, grouped by common prefixes """

      # Split parameter paths at '.' character
      paths = [ p.split('.', 1) for p in paths ]

      # Identify all names that are directly on this level
      # or represent parameter groups with a common prefix
      top_elements = { "%s%s" % (x[0], "=" if len(x) == 1 else ".") for x in paths }

      # Partition all names that are further down the tree by their prefix
      subpaths = {}
      for p in paths:
        if len(p) > 1:
          if p[0] not in subpaths: subpaths[p[0]] = []
          subpaths[p[0]].append(p[1])

      # If there are prefixes with only one name beneath them, put them on the top level
      for s in list(subpaths.iterkeys()):
        if len(subpaths[s]) == 1:
          top_elements.remove("%s." % s)
          top_elements.add("%s.%s=" % (s, subpaths[s][0]))
          del subpaths[s]

      result = { '': top_elements }
      # Revursively process each group
      for n, x in subpaths.iteritems():
        result[n] = construct_completion_tree(x)

      return result

    print 'function _dials_autocomplete_flags ()'
    print '{'
    print ' case "$1" in'
    for p in parameter_choice_list.iterkeys():
      print '\n  %s)' % p
      print '   _dials_autocomplete_values="%s";;' % ' '.join(parameter_choice_list[p])
    print '\n  *)'
    print '    _dials_autocomplete_values="";;'
    print ' esac'
    print '}'

    print 'function _dials_autocomplete_expansion ()'
    print '{'
    print ' case "$1" in'
    for p, exp in parameter_expansion_list.iteritems():
      if exp is not None:
        print '\n  %s=)' % p
        print '   _dials_autocomplete_values="%s=";;' % exp
    print '\n  *)'
    print '    _dials_autocomplete_values="";;'
    print ' esac'
    print '}'

    tree = construct_completion_tree(parameter_list)

    def _tree_to_bash(prefix, tree):
      for subkey in tree.iterkeys():
        if subkey != '':
          _tree_to_bash(prefix + subkey + '.', tree[subkey])
          print '\n  %s*)' % (prefix + subkey + '.')
          print '    _dials_autocomplete_values="%s";;' % " ".join(sorted([prefix + subkey + '.' + x for x in tree[subkey]['']]))

    print 'function _dials_autocomplete_hints ()'
    print '{'
    print ' case "$1" in'
    _tree_to_bash('', tree)

    toplevelset = tree[''] | set([p + "=" for p, exp in parameter_expansion_list.iteritems() if exp is not None])

    print '\n  *)'
    print '    _dials_autocomplete_values="%s";;' % " ".join(sorted(toplevelset))
    print ' esac'
    print '}'


def flatten_reflections(filename_object_list):
  '''
  Flatten a list of reflections tables

  :param filename_object_list: The parameter item
  :return: The flattened reflection table

  '''
  result = []
  for i in range(len(filename_object_list)):
    result.append(filename_object_list[i].data)
  return result

def flatten_datablocks(filename_object_list):
  '''
  Flatten a list of datablocks

  :param filename_object_list: The parameter item
  :return: The flattened list of datablocks

  '''
  result = []
  for i in range(len(filename_object_list)):
    result.extend(filename_object_list[i].data)
  return result

def flatten_experiments(filename_object_list):
  '''
  Flatten a list of experiment lists

  :param filename_object_list: The parameter item
  :return: The flattened experiment lists

  '''
  from dxtbx.model.experiment_list import ExperimentList
  result = ExperimentList()
  for i in range(len(filename_object_list)):
    result.extend(filename_object_list[i].data)
  return result
