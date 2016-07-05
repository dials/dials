#!/usr/bin/env python
#
# cluster.dials.exec
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
# LIBTBX_SET_DISPATCHER_NAME cluster.dials.exec

from __future__ import division

def get_cwd():
  '''
  Get the current working directory

  '''
  import sys
  return sys.argv[1]

def get_tid():
  '''
  Get the task id

  '''
  import os
  if "SGE_TASK_ID" in os.environ:
    return os.environ['SGE_TASK_ID']
  else:
    raise KeyError('Could not find task id')


if __name__ == '__main__':
  import traceback
  from os.path import join, exists
  from time import sleep
  import cPickle as pickle

  # Get the task id and the current working directory
  tid = get_tid()
  cwd = get_cwd()

  # Set the paths
  input_fn = join(cwd, "%s.input" % tid)
  output_fn = join(cwd, "%s.output" % tid)

  # Wait until it exists
  while not exists(input_fn):
    sleep(1)

  # Unpickle the function and the input
  with open(input_fn, "rb") as infile:
    function, element = pickle.load(open(input_fn))

  # Try to run the function, otherwise return an exception
  try:
    result = function(element)
  except Exception, exception:
    exception.args = [traceback.format_exc()]
    result = exception

  # Dump the result
  with open(output_fn, "wb") as outfile:
    pickle.dump(result, outfile, protocol=pickle.HIGHEST_PROTOCOL)
