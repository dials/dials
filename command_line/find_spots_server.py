from __future__ import division
import time
import BaseHTTPServer as server_base
from multiprocessing import Process as new_process
from multiprocessing import current_process
import os

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
  from dials.algorithms.peak_finding import per_image_analysis
  imageset = datablock.extract_imagesets()[0]
  stats = per_image_analysis.stats_single_image(imageset, reflections)
  return stats
  return stats.n_spots_total, stats.n_spots_no_ice, #stats.n_spots_4A, stats.total_intensity, stats.estimated_d_min

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
    proc = current_process().name
    try:
      stats = work(filename, params)
      s.wfile.write("<response>%s: %6d / %6d</response>" %
                    (os.path.split(filename)[-1],
                     stats.n_spots_total, stats.n_spots_no_ice))

    except:
      s.wfile.write("<response>error</response>")
    return

def serve(httpd):
  try:
    while not stop:
      httpd.handle_request()
  except KeyboardInterrupt:
    pass
  return

def main(nproc):
  server_class = server_base.HTTPServer
  httpd = server_class(('', 1701), handler)
  print time.asctime(), 'start'
  for j in range(nproc - 1):
    new_process(target=serve, args=(httpd,)).start()
  serve(httpd)
  httpd.server_close()
  print time.asctime(), 'done'

if __name__ == '__main__':
  import sys
  if len(sys.argv) > 1:
    nproc = int(sys.argv[1])
  else:
    from libtbx.introspection import number_of_processors
    nproc = number_of_processors(return_value_if_unknown=-1)
  main(nproc)
