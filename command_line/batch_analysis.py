from __future__ import division

import os

from cctbx.array_family import flex
import iotbx.phil


phil_scope = iotbx.phil.parse('''\
nproc = Auto
  .type = int(value_min=1)
json = None
  .type = path
''')

help_message = '''\
'''

def work(args):
  filename = args[0]
  cl = args[1]
  from dials.command_line import find_spots_server

  d = {'image': filename}
  d.update(find_spots_server.work(filename, cl=cl))

  return d

def work_all(filenames, args, nproc):
  from libtbx import easy_mp

  cl = args
  args = []
  for f in filenames:
    args.append((f, cl))

  results = easy_mp.parallel_map(
    func=work,
    iterable=args,
    processes=nproc,
    method="multiprocessing",
    preserve_order=True,
    asynchronous=True,
    preserve_exception_message=True)

  return results


def run(args):
  from dials.util.options import OptionParser
  import libtbx.load_env

  usage = "%s [options] find_spots.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    epilog=help_message)

  filenames = [arg for arg in args if os.path.isfile(arg)]
  args = [arg for arg in args if not arg in filenames]

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)
  if params.nproc is libtbx.Auto:
    from libtbx.introspection import number_of_processors
    params.nproc = number_of_processors(return_value_if_unknown=-1)
  print 'nproc: %i' %params.nproc
  results = work_all(filenames, args, nproc=params.nproc)
  print results

  if params.json is not None:
    import json
    with open(params.json, 'wb') as f:
      json.dump(results, f)


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
