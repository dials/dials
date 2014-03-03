#!/usr/bin/env python
#
# logging.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from libtbx.phil import parse

phil_scope = parse('''

logging {
  version = 1
    .type = int
  disable_existing_loggers = True
    .type = bool

  formatters {
    standard {
      format = '%(asctime)s [%(levelname)s] %(name)s: %(message)s'
        .type = str
    }
  }

  handlers {
    console {
      level = DEBUG INFO *WARN ERROR CRITICAL NOTSET
        .type = choice
      class = 'logging.StreamHandler'
        .type = str
      formatter = *standard
        .type = choice
    }
    file {
      level = DEBUG *INFO WARN ERROR CRITICAL NOTSET
        .type = choice
      class = 'logging.FileHandler'
        .type = str
      formatter = *standard
        .type = choice
      directory = '~/.dials'
        .type = str
    }
  }

  loggers {
    default {
      handlers = console file
        .type = choice(multi=True)
      level = *DEBUG INFO WARN ERROR CRITICAL NOTSET
        .type = choice
      propagate = True
        .type = bool
    }
  }
}

''')
