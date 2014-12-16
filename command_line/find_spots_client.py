from __future__ import division

def work(host, filename, params):
  import httplib
  conn = httplib.HTTPConnection(host, 1701)
  path = filename
  for param in params:
    path += ';%s' % param
  conn.request('GET', path)
  response = conn.getresponse()
  return response.read()

def nproc():
  from libtbx.introspection import number_of_processors
  return number_of_processors(return_value_if_unknown=-1)

def work_all(host, filenames, params):
  from multiprocessing.pool import ThreadPool as thread_pool
  pool = thread_pool(processes=nproc())
  threads = { }
  for filename in filenames:
    threads[filename] = pool.apply_async(work, (host, filename, params))
  results = { }
  for filename in filenames:
    results[filename] = threads[filename].get()
    print results[filename]
  return

def stop(host, nproc):
  import httplib
  from socket import error as socket_error
  for j in range(nproc):
    try:
      conn = httplib.HTTPConnection(host, 1701)
      path = '/Ctrl-C'
      conn.request('GET', path)
      response = conn.getresponse()
    except socket_error, e:
      # run out of procs
      break
  return j

if __name__ == '__main__':
  import sys

  if len(sys.argv) < 3:
    raise RuntimeError, '%s [host] [filename] [param=value]'

  if sys.argv[2] == 'stop':
    if len(sys.argv) < 3:
      raise RuntimeError, '%s stop [nproc]'
    host = sys.argv[1]
    if len(sys.argv) > 3:
      nproc = int(sys.argv[3])
    else:
      nproc = 1024
    stopped = stop(host, nproc)
    print 'Stopped %d findspots processes' % stopped
  else:
    host = sys.argv[1]
    filenames = []
    params = []
    for arg in sys.argv[2:]:
      if '=' in arg:
        params.append(arg)
      else:
        filenames.append(arg)
    if len(filenames) == 1:
      print work(host, filenames[0], params)
    else:
      work_all(host, filenames, params)
