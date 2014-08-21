#!/usr/bin/env python
#
# indexing.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

def generate_phil_scope():

  from libtbx.phil import parse
  from dials.algorithms.indexing.indexer import master_phil_scope

  phil_scope = parse('''

  indexing
    .help = "Configure the indexing algorithm."
  {
  }

  ''', process_includes=True)

  main_scope = phil_scope.get_without_substitution("indexing")
  assert(len(main_scope) == 1)
  main_scope = main_scope[0]
  main_scope.adopt_scope(master_phil_scope)

  return phil_scope

phil_scope = generate_phil_scope()

