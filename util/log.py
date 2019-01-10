#!/usr/bin/env python
#
# log.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import logging


def config(verbosity=1, name='dials', info=None, debug=None):
  """
  Configure the logging.

  :param verbosity: Verbosity level of stdout log output.  Possible values:
                      * 0: No log output to stdout;
                      * 1: Info log output to stdout;
                      * 2: Info & debug log output to stdout.
  :type verbosity: int
  :param name: Logger name.
  :type name: str
  :param info: Filename for info log output.  If False, no info log file is
               written.
  :type info: str
  :param debug: Filename for debug log output.  If False, no debug log file is
                written.
  :type debug: str
  """

  import logging.config

  # Debug or not
  if verbosity > 1:
    level = 'DEBUG'
  else:
    level = 'INFO'

  # Preapare a dictionary of FileHandler configuration options
  handlers_dict = {
    'stream' : {
      'level' : level,
      'class' : 'logging.StreamHandler',
      'formatter' : 'standard',
      'stream': 'ext://sys.stdout',
    }
  }
  # Set the handlers to use
  if verbosity > 0:
    handlers = ['stream']
  else:
    handlers = []
  if info:
    # Prepare a FileHandler for the info log
    handlers_dict.update({
      'file_info' : {
        'level' : 'INFO',
        'class' : 'logging.FileHandler',
        'formatter' : 'standard',
        'filename' : info,
        'mode': 'w'
      }
    })
    # Add the info log FileHandler to the list in the Logger
    handlers.append('file_info')
  if debug:
    # Prepare a FileHandler for the debug log
    handlers_dict.update({
      'file_debug' : {
        'level' : 'DEBUG',
        'class' : 'logging.FileHandler',
        'formatter' : 'standard',
        'filename' : debug,
        'mode': 'w'
      },
    })
    # Add the debug log FileHandler to the Logger
    handlers.append('file_debug')

  # Configure the logging
  config_dict = {
    'version' : 1,
    'disable_existing_loggers' : False,

    'formatters' : {
      'standard' : {
        'format' : '%(message)s'
      },
      'extended' : {
        'format' : '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
      },
    },

    'handlers' : handlers_dict,

    'loggers' : {
      name : {
        'handlers' : handlers,
        'level' : 'DEBUG',
        'propagate' : True
      },
    }
  }
  logging.config.dictConfig(config_dict)
  logging.getLogger('dxtbx').setLevel(logging.DEBUG)
  import dials.util.banner


def config_simple_stdout(name='dials'):
  '''
  Configure the logging to just go to stdout

  '''
  import logging.config

  # Configure the logging
  logging.config.dictConfig({

    'version' : 1,
    'disable_existing_loggers' : False,

    'formatters' : {
      'standard' : {
        'format' : '%(message)s'
      }
    },

    'handlers' : {
      'stream' : {
        'level' : 'DEBUG',
        'class' : 'logging.StreamHandler',
        'formatter' : 'standard',
        'stream': 'ext://sys.stdout',
      }
    },

    'loggers' : {
      name : {
        'handlers' : ['stream'],
        'level' : 'DEBUG',
        'propagate' : True
      },
    }
  })

  import dials.util.banner


class CacheHandler(logging.Handler):
  ''' A simple class to store log messages. '''

  def __init__(self):
    '''
    Initialise the handler

    '''
    super(CacheHandler, self).__init__()
    self._messages = []

  def emit(self, record):
    '''
    Emit the message to a list

    :param record: The log record

    '''
    self._messages.append(record)

  def messages(self):
    return self._messages


def config_simple_cached():
  '''
  Configure the logging to use a cache.

  '''
  import logging.config

  # Configure the logging
  logging.config.dictConfig({

    'version' : 1,
    'disable_existing_loggers' : False,

    'handlers' : {
      'cache' : {
        'level' : 'DEBUG',
        'class' : 'dials.util.log.CacheHandler',
      }
    },

    'loggers' : {
      'dials' : {
        'handlers' : ['cache'],
        'level' : 'DEBUG',
        'propagate' : True
      }
    }
  })


class LoggerIO(object):
  ''' Wrap the logger with file type object '''

  def __init__(self, logger, level):
    '''
    Initialise the logger io

    :param level: The logging level

    '''
    self.logger = logger
    self.level = level

  def write(self, buf):
    '''
    Write to the logger

    :param buf: The buffer

    '''
    self.logger.log(self.level, buf)

  def flush(self):
    '''
    Flush (don't do anything)

    '''
    pass


def info_handle(logger):
  '''
  :return: A handle to an INFO logger file object

  '''
  return LoggerIO(logger, logging.INFO)


def debug_handle(logger):
  '''
  :return: A handle to an DEBUG logger file object

  '''
  return LoggerIO(logger, logging.DEBUG)
