from __future__ import division
import time
import BaseHTTPServer as server_base

stop = False

def work(filename, cl=[]):
  from dials.command_line.find_spots import phil_scope as params
  from dxtbx.datablock import DataBlockFactory
  from dials.array_family import flex
  interp = params.command_line_argument_interpreter()
  for cla in cl:
    params = params.fetch(interp.process(cla))
  datablock = DataBlockFactory.from_filenames([filename])[0]
  reflections = flex.reflection_table.from_observations(
    datablock, params.extract())
  return len(reflections)

class handler(server_base.BaseHTTPRequestHandler):
  def do_GET(s):
    """Respond to a GET request."""
    s.send_response(200)
    s.send_header("Content-type", "text/xml")
    s.end_headers()
    if s.path == '/Ctrl-C':
      global stop
      stop = True
      return
    filename = s.path.split(';')[0]
    params = s.path.split(';')[1:]
    try:
      result = work(filename, params)
      s.wfile.write("<response>%s: %d</response>" % (filename, result))
    except:
      s.wfile.write("<response>error</response>")
    return

if __name__ == '__main__':
  server_class = server_base.HTTPServer
  httpd = server_class(('', 1701), handler)
  print time.asctime(), 'start'
  try:
    while not stop:
      httpd.handle_request()
  except KeyboardInterrupt:
    pass
  httpd.server_close()
  print time.asctime(), 'done'
