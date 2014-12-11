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

from __future__ import division

def config(verbosity=1, filename=''):
  ''' Configure the logging. '''
  import logging.config

  # Debug or not
  if verbosity > 1:
    logging_level = 'DEBUG'
  else:
    logging_level = 'INFO'

  # Set the handlers to use
  if verbosity > 0:
    handlers = ['stream']
  else:
    handlers = []
  if filename is not None and filename != '':
    handlers.append('file')
  else:
    filename = 'dials.log'

  # Configure the logging
  logging.config.dictConfig({

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

    'handlers' : {
      'stream' : {
        'level' : 'DEBUG',
        'class' : 'logging.StreamHandler',
        'formatter' : 'standard',
      },
      'file' : {
        'level' : 'DEBUG',
        'class' : 'logging.FileHandler',
        'formatter' : 'standard',
        'filename' : filename,
        'mode' : 'w'
      }
    },

    'loggers' : {
      '' : {
        'handlers' : handlers,
        'level' : logging_level,
        'propagate' : True
      }
    }
  })

def config_simple_stdout():
  ''' Configure the logging. '''
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
      }
    },

    'loggers' : {
      '' : {
        'handlers' : ['stream'],
        'level' : 'DEBUG',
        'propagate' : True
      }
    }
  })

