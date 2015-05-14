from __future__ import division
import time
import BaseHTTPServer as server_base
from multiprocessing import Process as new_process
from multiprocessing import current_process
import os

stop = False

def work(filename, cl=[]):
  import libtbx.phil
  phil_scope = libtbx.phil.parse('''\
index = False
  .type = bool
''')
  interp = phil_scope.command_line_argument_interpreter()
  params, unhandled = interp.process_and_fetch(
    cl, custom_processor='collect_remaining')
  index = params.extract().index

  from dials.command_line.find_spots import phil_scope as params
  from dxtbx.datablock import DataBlockFactory
  from dials.array_family import flex
  interp = params.command_line_argument_interpreter()
  params, unhandled = interp.process_and_fetch(
    unhandled, custom_processor='collect_remaining')
  datablock = DataBlockFactory.from_filenames([filename])[0]
  reflections = flex.reflection_table.from_observations(
    datablock, params.extract())
  from dials.algorithms.peak_finding import per_image_analysis
  imageset = datablock.extract_imagesets()[0]
  stats = per_image_analysis.stats_single_image(
    imageset, reflections,
    i=imageset.get_scan().get_image_range()[0]-1, plot=False)

  if index and stats.n_spots_no_ice > 10:
    import logging
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    from dials.command_line.index import phil_scope
    interp = phil_scope.command_line_argument_interpreter()
    params, unhandled = interp.process_and_fetch(
      unhandled, custom_processor='collect_remaining')
    imagesets = [imageset]
    params = params.extract()
    params.indexing.scan_range=[]
  
    if (imageset.get_goniometer() is not None and
        imageset.get_scan() is not None and
        imageset.get_scan().get_oscillation()[1] == 0):
      imageset.set_goniometer(None)
      imageset.set_scan(None)

    try:
      if params.indexing.method == "fft3d":
        from dials.algorithms.indexing.fft3d import indexer_fft3d
        idxr = indexer_fft3d(reflections, imagesets, params=params)
      elif params.indexing.method == "fft1d":
        from dials.algorithms.indexing.fft1d import indexer_fft1d
        idxr = indexer_fft1d(reflections, imagesets, params=params)
      elif params.indexing.method == "real_space_grid_search":
        from dials.algorithms.indexing.real_space_grid_search \
             import indexer_real_space_grid_search
        idxr = indexer_real_space_grid_search(reflections, imagesets, params=params)
      stats.crystal = idxr.refined_experiments.crystals()[0]
      stats.n_indexed = len(idxr.refined_reflections)
    except Exception, e:
      stats.crystal = None
      stats.n_indexed = None
    print logging

  return stats
  return stats.n_spots_total, stats.n_spots_no_ice

class handler(server_base.BaseHTTPRequestHandler):
  def do_GET(s):
    '''Respond to a GET request.'''
    s.send_response(200)
    s.send_header('Content-type', 'text/xml')
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
      response = [
        '<image>%s</image>' % filename,
        '<spot_count>%s</spot_count>' % stats.n_spots_total,
        '<spot_count_no_ice>%s</spot_count_no_ice>' % stats.n_spots_no_ice,
        '<d_min>%.2f</d_min>' % stats.estimated_d_min,
        '<d_min_method_1>%.2f</d_min_method_1>' % stats.d_min_distl_method_1,
        '<d_min_method_2>%.2f</d_min_method_2>' % stats.d_min_distl_method_2,
        '<total_intensity>%.0f</total_intensity>' % stats.total_intensity,
      ]
      if hasattr(stats, 'crystal') and stats.crystal is not None:
        response.append(
          '<unit_cell>%.6g %.6g %.6g %.6g %.6g %.6g</unit_cell>' %stats.crystal.get_unit_cell().parameters())
        response.append(
          '<n_indexed>%i</n_indexed>' %stats.n_indexed)
        
      s.wfile.write('<response>\n%s\n</response>' % ('\n'.join(response)))

    except Exception, e:
      #import traceback
      #traceback.print_exc()
      s.wfile.write('<response>error: %s</response>' % str(e))
    return

def serve(httpd):
  try:
    while not stop:
      httpd.handle_request()
  except KeyboardInterrupt:
    pass
  return


import libtbx.phil
phil_scope = libtbx.phil.parse('''\
nproc = Auto
  .type = int(value_min=2)
port = 1701
  .type = int(value_min=1)
''')


def main(nproc, port):
  server_class = server_base.HTTPServer
  httpd = server_class(('', port), handler)
  print time.asctime(), 'start'
  for j in range(nproc - 1):
    new_process(target=serve, args=(httpd,)).start()
  serve(httpd)
  httpd.server_close()
  print time.asctime(), 'done'

if __name__ == '__main__':
  import sys
  import libtbx.load_env

  usage = '%s [options]' % libtbx.env.dispatcher_name

  from dials.util.options import OptionParser
  parser = OptionParser(usage=usage, phil=phil_scope)
  params, options = parser.parse_args(show_diff_phil=True)
  if params.nproc is libtbx.Auto:
    from libtbx.introspection import number_of_processors
    params.nproc = number_of_processors(return_value_if_unknown=-1)
  main(params.nproc, params.port)
