import httplib

def work(filename, params):
  conn = httplib.HTTPConnection('localhost', 1701)
  path = filename
  for param in params:
    path += ';%s' % param
  conn.request('GET', path)
  response = conn.getresponse()
  return response.read()

if __name__ == '__main__':
  import sys
  filename = sys.argv[1]
  params = sys.argv[2:]
  print work(filename, params)
