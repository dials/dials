#
# command_line.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import division

def parse_range_list_string(string):
  """Parse a string in the following ways:
  string: 1, 2, 3        -> [1, 2, 3]
  string: 1 - 6          -> [1, 2, 3, 4, 5, 6]
  string: 1 - 6, 7, 8, 9 -> [1, 2, 3, 4, 5, 6, 7, 8, 9]
  """
  items = string.split(',')
  for i in range(len(items)):
    items[i] = items[i].split("-")
    if len(items[i]) == 1:
      items[i] = [int(items[i][0])]
    elif len(items[i]) == 2:
      items[i] = range(int(items[i][0]), int(items[i][1]) + 1)
    else:
      raise SyntaxError
  items = [item for sublist in items for item in sublist]
  return set(items)

def interactive_console():
  """ Enter an interactive console session. """
  try:
    from IPython import embed
    import inspect
    frame = inspect.currentframe()
    try:
      embed(user_ns = frame.f_back.f_locals)
    finally:
        del frame
  except ImportError:
    print "IPython not available"


class ProgressBarTimer:
  """ A simple timer for the progress bar. """

  def __init__(self):
    """ Init the progress bar timer. """
    from time import time
    self._last_time = time()
    self._last_perc = 0
    self._update_period = 0.5
    self._n_seconds_left = -1

  def update(self, percent):
    """ Update the timer. """
    from time import time

    # Get the current time diff between last time
    curr_time = time()
    diff_time = curr_time - self._last_time

    # Only update after certain period or at 100%
    if percent < 0: percent = 0
    if percent > 100: percent = 100
    if diff_time >= self._update_period or percent >= 100:

      # Check the difference in percentage and calculate
      # number of seconds remaining
      diff_perc = percent - self._last_perc
      if (diff_perc == 0):
        self._n_seconds_left = 0
      else:
        self._n_seconds_left = diff_time * (100 - percent) / diff_perc

    # Return number of seconds
    return self._n_seconds_left

class ProgressBar:
  """ A command line progress bar. """

  def __init__(self, title=None, spinner=True, bar=True, estimate_time=True,
               indent=0, length=80):
    """ Init the progress bar parameters. """

    # Set the parameters
    self._title = title
    self._indent = indent
    self._spinner = spinner
    self._estimate_time = estimate_time
    self._bar = bar
    self._length = length

    self._timer = ProgressBarTimer()
    self._start_time = self._timer._last_time

    # Print 0 percent
    self.update(0)

  def update(self, fpercent):
    """ Update the progress bar with a percentage. """
    import sys
    from math import ceil

    # do not update if not a tty
    if not sys.stdout.isatty():
      return

    # Get integer percentage
    percent = int(fpercent)
    if percent < 0: percent = 0
    if percent > 100: percent = 100

    # Add a percentage counter
    right_str = ''
    left_str = ''
    if sys.stdout.isatty():
      left_str = '\r'
    left_str += ' ' * self._indent

    # Add a title if given
    if self._title:
      left_str += self._title + ': '

    left_str += '{0: >3}%'.format(percent)

    # Add a spinner
    if self._spinner:
      left_str += ' '
      left_str += '[ {0} ]'.format('-\|/'[percent % 4])

    # Add a timer
    if self._estimate_time:
      n_seconds_left = self._timer.update(fpercent)
      if n_seconds_left < 0:
        n_seconds_left = '?'
      else:
        n_seconds_left = int(ceil(n_seconds_left))
      right_str = ' ' + 'est: {0}s'.format(n_seconds_left) + right_str

    # Add a bar
    if self._bar:
      bar_length = self._length - (len(left_str) + len(right_str)) - 5
      n_char = int(percent * bar_length / 100)
      n_space = bar_length - n_char
      left_str += ' '
      left_str += '[ {0}>{1} ]'.format('=' * n_char, ' ' * n_space)

    # Append strings
    progress_str = left_str + right_str

    # Print progress string to stdout
    sys.stdout.write(progress_str)
    sys.stdout.flush()

  def finished(self, string=None):
    """ The progress bar is finished. """
    if string:
      self._title = string
    else:
      string = ''

    ''' Print the 'end of comand' string.'''
    from sys import stdout
    from time import time


    if self._estimate_time:
      # Get the time string
      time_string = '{0:.2f}s'.format(time() - self._start_time)

      # Truncate the string
      max_length = self._length - self._indent - len(time_string) - 1
      string = string[:max_length]

      # Add an indent and a load of dots and then the time string
      dot_length = 1 + max_length - len(string)
      string = (' ' * self._indent) + string
      string = string + '.' * (dot_length)
      string = string + time_string

    else:

      # Truncaet the string
      max_length = self._length - self._indent
      string = string[:max_length]

      # Add a load of dots
      dot_length = max_length - len(string)
      string = (' ' * self._indent) + string
      string = string + '.' * (dot_length)

    # Write the string to stdout
    if stdout.isatty():
      string = '\r' + string + '\n'
    else:
      string = string + '\n'
    stdout.write(string)
    stdout.flush()

class Command(object):
  '''Class to nicely print out a command with timing info.'''

  # Variables available in class methods
  indent = 0
  max_length = 80
  print_time = True

  @classmethod
  def start(self, string):
    ''' Print the 'start command' string.'''
    from sys import stdout
    from time import time
    # from termcolor import colored

    # Get the command start time
    self._start_time = time()

    # do not output if not a tty
    if not stdout.isatty():
      return

    # Truncate the string to the maximum length
    max_length = self.max_length - self.indent - 3
    string = string[:max_length]
    string = (' ' * self.indent) + string + '...'

    # Write the string to stdout
    stdout.write(string)
    stdout.flush()

  @classmethod
  def end(self, string):
    ''' Print the 'end of comand' string.'''
    from sys import stdout
    from time import time
    #from termcolor import colored

    # Check if we want to print the time or not
    if self.print_time:

      # Get the time string
      time_string = '{0:.2f}s'.format(time() - self._start_time)

      # Truncate the string
      max_length = self.max_length - self.indent - len(time_string) - 1
      string = string[:max_length]

      # Add an indent and a load of dots and then the time string
      dot_length = 1 + max_length - len(string)
      string = (' ' * self.indent) + string
      string = string + '.' * (dot_length)
      string = string + time_string

    else:

      # Truncaet the string
      max_length = self.max_length - self.indent
      string = string[:max_length]

      # Add a load of dots
      dot_length = max_length - len(string)
      string = (' ' * self.indent) + string
      string = string + '.' * (dot_length)

    # Write the string to stdout
    if stdout.isatty():
      string = '\r' + string + '\n'
    else:
      string = string + '\n'
    stdout.write(string)
    stdout.flush()


class Importer(object):
  ''' A class to import the command line arguments. '''

  def __init__(self, args, include=None, exclude=None, verbose=False,
               check_format=True):
    ''' Parse the arguments.

    If include is set, only those items set will be tried. If not, then if
    exclude is set, then those items will not be tested.

    These are the types we can import:
     - images: a list of images
     - reflections : a list of reflections
     - datablocks : a list of datablocks
     - experiments: a list of experiments

    Params:
      args The arguments to parse
      include types to try
      exclude types not to try
      verbose True/False print out some stuff

    Example:
      import = Importer(argv, include=['reflections'])

    '''

    # Check the format in data block and experiment list
    self._check_format = check_format

    # Initialise output
    self.datablocks = None
    self.experiments = None
    self.reflections = None

    # Get the list of items to try
    totry = ['pickle', 'images', 'datablocks',
             'experiments', 'reflections']
    if include is not None:
      for item in include:
        assert(item in totry)
      totry = include
    elif exclude is not None:
      for item in exclude:
        totry.remove(item)

    # Try to import the items
    unhandled = args
    for item in totry:
      if verbose: print 'Try import as %s' % item
      unhandled = self.try_import(unhandled, item, verbose)
    self.unhandled_arguments = unhandled

  def try_import_pickle(self, args, verbose, pickle_ext=('.pickle', '.pkl')):
    import os
    from libtbx import easy_pickle
    unhandled = []
    for i, arg in enumerate(args):
      if isinstance(arg, basestring) and os.path.splitext(arg)[1] in pickle_ext:
        try:
          from dials.array_family import flex
          obj = easy_pickle.load(arg)
          unhandled.append(obj)
        except Exception:
          unhandled.append(arg)
      else:
        unhandled.append(arg)
    return unhandled

  def try_import(self, args, item, verbose):
    ''' Try to import with the given item. '''
    return getattr(self, "try_import_%s" % item)(args, verbose)

  def try_import_images(self, args, verbose):
    ''' Try to import images. '''
    from dxtbx.datablock import DataBlockFactory
    unhandled = []
    datablocks = DataBlockFactory.from_args(
      args, verbose=verbose, unhandled=unhandled)
    if len(datablocks) > 0:
      if self.datablocks == None:
        self.datablocks = datablocks
      else:
        self.datablocks.extend(datablocks)
    return unhandled

  def try_import_datablocks(self, args, verbose):
    ''' Try to import imagesets. '''
    from dxtbx.datablock import DataBlockFactory, DataBlock
    from os.path import abspath
    unhandled = []
    for argument in args:
      if isinstance(argument, DataBlock):
        if self.datablocks == None:
          self.datablocks = [argument]
        else:
          self.datablocks.append(argument)
        continue
      try:
        abs_argument = abspath(argument)
        datablocks = DataBlockFactory.from_serialized_format(
          abs_argument, self._check_format)
        if verbose: print 'Loaded %s as datablock list' % abs_argument
        if self.datablocks == None:
          self.datablocks = datablocks
        else:
          self.datablocks.extend(datablocks)
      except Exception:
        unhandled.append(argument)
    return unhandled

  def try_import_experiments(self, args, verbose):
    ''' Try to import experiments. '''
    from dxtbx.model.experiment.experiment_list import ExperimentListFactory
    from os.path import abspath
    unhandled = []
    for argument in args:
      try:
        abs_argument = abspath(argument)
        experiments = ExperimentListFactory.from_serialized_format(
          abs_argument, self._check_format)
        if verbose: print 'Loaded %s as experiment list' % abs_argument
        if self.experiments == None:
          self.experiments = experiments
        else:
          self.experiments.extend(experiments)
      except Exception:
        unhandled.append(argument)
    return unhandled

  def try_import_reflections(self, args, verbose):
    ''' Try to import reflections. '''
    from dials.array_family import flex
    import cPickle as pickle
    unhandled = []
    for argument in args:
      try:
        if isinstance(argument, flex.reflection_table):
          obj = argument
        else:
          with open(argument, 'rb') as inputfile:
            obj = pickle.load(inputfile)
            if not isinstance(obj, flex.reflection_table):
              continue
        if verbose:
          print 'Loaded %s as reflection table' % argument
        if self.reflections is None:
          self.reflections = [obj]
        else:
          self.reflections.append(obj)
      except Exception:
        unhandled.append(argument)
    return unhandled


if __name__ == '__main__':
  import time

  p = ProgressBar()

  for j in range(100):
    p.update(j)
    time.sleep(0.05)

  p.finished()

  Command.start("Starting to do a command")
  time.sleep(1)
  Command.end("Ending the command")
