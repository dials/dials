from __future__ import division
import time
import BaseHTTPServer as server_base

stop = False

class handler(server_base.BaseHTTPRequestHandler):
  def do_GET(s):
    """Respond to a GET request."""
    s.send_response(200)
    s.send_header("Content-type", "text/xml")
    s.end_headers()
    if s.path == '/Ctrl-C':
      global stop
      stop = True
    s.wfile.write("<response>path: %s</response>" % s.path)
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
