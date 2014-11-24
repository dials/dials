import httplib

def work(filename, params):
  conn = httplib.HTTPConnection('localhost', 1701)
  path = filename
  for param in params:
    path += ';%s' % param
  conn.request('GET', path)
  response = conn.getresponse()
  return response.read()

def stop(nproc):
  for j in range(nproc):
    conn = httplib.HTTPConnection('localhost', 1701)
    path = '/Ctrl-C'
    conn.request('GET', path)
    response = conn.getresponse()
  return


if __name__ == '__main__':
  import sys

  if len(sys.argv) < 2:
    raise RuntimeError, '%s [filename] [param=value]'

  if sys.argv[1] == 'stop':
    if len(sys.argv) < 3:
      raise RuntimeError, '%s stop nproc'
    nproc = int(sys.argv[2])
    stop(nproc)
  else:
    filename = sys.argv[1]
    params = sys.argv[2:]
    print work(filename, params)
