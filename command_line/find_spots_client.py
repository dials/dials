from __future__ import division

def work(host, port, filename, params):
  import httplib
  conn = httplib.HTTPConnection(host, port)
  path = filename
  for param in params:
    path += ';%s' % param
  conn.request('GET', path)
  response = conn.getresponse()
  return response.read()

def nproc():
  from libtbx.introspection import number_of_processors
  return number_of_processors(return_value_if_unknown=-1)

def work_all(host, port, filenames, params):
  from multiprocessing.pool import ThreadPool as thread_pool
  pool = thread_pool(processes=nproc())
  threads = { }
  for filename in filenames:
    threads[filename] = pool.apply_async(work, (host, port, filename, params))
  results = { }
  for filename in filenames:
    results[filename] = threads[filename].get()
    print results[filename]
  return

def stop(host, port, nproc):
  import httplib
  from socket import error as socket_error
  for j in range(nproc):
    try:
      conn = httplib.HTTPConnection(host, port)
      path = '/Ctrl-C'
      conn.request('GET', path)
      response = conn.getresponse()
    except socket_error, e:
      # run out of procs
      break
  return j

import libtbx.phil
phil_scope = libtbx.phil.parse("""\
nproc = Auto
  .type = int(value_min=2)
host = localhost
  .type = str
port = 1701
  .type = int(value_min=1)
""")

if __name__ == '__main__':
  import os
  import sys
  import libtbx.load_env

  args = sys.argv[1:]

  filenames = [arg for arg in args if os.path.isfile(arg)]
  args = [arg for arg in args if not arg in filenames]

  usage = "%s [options] filenames" %libtbx.env.dispatcher_name

  from dials.util.options import OptionParser
  parser = OptionParser(usage=usage, phil=phil_scope)
  params, options, args = parser.parse_args(
    args=args, show_diff_phil=True, return_unhandled=True)
  if params.nproc is libtbx.Auto:
    params.nproc = 1024

  if args[0] == 'stop':
    stopped = stop(params.host, params.port, params.nproc)
    print 'Stopped %d findspots processes' % stopped
  else:
    if len(filenames) == 1:
      print work(params.host, params.port, filenames[0], args)
    else:
      work_all(params.host, params.port, filenames, args)
