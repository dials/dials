import httplib

def work(host, filename, params):
  conn = httplib.HTTPConnection(host, 1701)
  path = filename
  for param in params:
    path += ';%s' % param
  conn.request('GET', path)
  response = conn.getresponse()
  return response.read()

def stop(host, nproc):
  for j in range(nproc):
    conn = httplib.HTTPConnection(host, 1701)
    path = '/Ctrl-C'
    conn.request('GET', path)
    response = conn.getresponse()
  return


if __name__ == '__main__':
  import sys

  if len(sys.argv) < 3:
    raise RuntimeError, '%s [host] [filename] [param=value]'

  if sys.argv[2] == 'stop':
    if len(sys.argv) < 4:
      raise RuntimeError, '%s stop nproc'
    host = sys.argv[1]
    nproc = int(sys.argv[3])
    stop(host, nproc)
  else:
    host = sys.argv[1]
    filename = sys.argv[2]
    params = sys.argv[3:]
    print work(host, filename, params)
