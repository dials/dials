from __future__ import division
import time
import BaseHTTPServer as server_base
from multiprocessing import Process
from multiprocessing import current_process
import os
import signal

help_message = '''\
A client/server version of dials.find_spots with additional analysis including
estimation of resolution limits. Intended for quick feedback of image quality
during grid scans and data collections.

On the server machine::

  dials.find_spots_server [nproc=8] [port=1234]

On the client machine::

  dials.find_spots_client [host=hostname] [port=1234] [nproc=8] /path/to/image.cbf

The client will return a short xml string indicating the number of spots found
and several estimates of the resolution limit.

e.g.::

  <response>
  <image>/path/to/image_0001.cbf</image>
  <spot_count>352</spot_count>
  <spot_count_no_ice>263</spot_count_no_ice>
  <d_min>1.46</d_min>
  <d_min_method_1>1.92</d_min_method_1>
  <d_min_method_2>1.68</d_min_method_2>
  <total_intensity>56215</total_intensity>
  </response>


* ``spot_count`` is the total number of spots found in given image
* ``spot_count_no_ice`` is the number of spots found excluding those at resolutions
  where ice rings may be found
* ``d_min_method_1`` is equivalent to distl's resolution estimate method 1
* ``d_min_method_2`` is equivalent to distl's resolution estimate method 2
* ``total_intensity`` is the total intensity of all strong spots excluding those
  at resolutions where ice rings may be found

Any valid ``dials.find_spots`` parameter may be passed to
``dials.find_spots_client``, e.g.::

  dials.find_spots_client /path/to/image.cbf min_spot_size=2 d_min=2

To stop the server::

  dials.find_spots_client stop [host=hostname] [port=1234]

'''

stop = False

def work(filename, cl=[]):
  import libtbx.phil
  phil_scope = libtbx.phil.parse('''\
index = False
  .type = bool
integrate = False
  .type = bool
''')
  interp = phil_scope.command_line_argument_interpreter()
  params, unhandled = interp.process_and_fetch(
    cl, custom_processor='collect_remaining')
  index = params.extract().index
  integrate = params.extract().integrate

  from dials.command_line.find_spots import phil_scope as params
  from dxtbx.datablock import DataBlockFactory
  from dials.array_family import flex
  interp = params.command_line_argument_interpreter()
  params, unhandled = interp.process_and_fetch(
    unhandled, custom_processor='collect_remaining')
  params = params.extract()
  # no need to write the hot mask in the server/client
  params.spotfinder.write_hot_mask = False
  datablock = DataBlockFactory.from_filenames([filename])[0]
  reflections = flex.reflection_table.from_observations(datablock, params)
  from dials.algorithms.peak_finding import per_image_analysis
  imageset = datablock.extract_imagesets()[0]
  stats = per_image_analysis.stats_single_image(
    imageset, reflections,
    i=imageset.get_scan().get_image_range()[0]-1, plot=False)

  if index and stats.n_spots_no_ice > 10:
    import logging
    logging.basicConfig(stream=sys.stdout, level=logging.INFO)
    from dials.algorithms.indexing import indexer
    interp = indexer.master_phil_scope.command_line_argument_interpreter()
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
      idxr = indexer.indexer_base.from_parameters(
        reflections, imagesets, params=params)
      stats.crystal = idxr.refined_experiments.crystals()[0]
      stats.n_indexed = len(idxr.refined_reflections)
      stats.fraction_indexed = stats.n_indexed/len(reflections)
    except Exception, e:
      print e
      stats.crystal = None
      stats.n_indexed = None
      stats.fraction_indexed = None

    if integrate and stats.crystal is not None:

      from dials.algorithms.profile_model.factory import ProfileModelFactory
      from dials.algorithms.integration.integrator import IntegratorFactory
      from dials.command_line.integrate import phil_scope
      interp = phil_scope.command_line_argument_interpreter()
      params, unhandled = interp.process_and_fetch(
        unhandled, custom_processor='collect_remaining')
      imagesets = [imageset]
      params = params.extract()

      params.profile.gaussian_rs.min_spots = 0

      experiments = idxr.refined_experiments
      reference = idxr.refined_reflections

      predicted = flex.reflection_table.from_predictions_multi(
        experiments,
        dmin=params.prediction.d_min,
        dmax=params.prediction.d_max,
        margin=params.prediction.margin,
        force_static=params.prediction.force_static
        )

      matched, reference = predicted.match_with_reference(reference)
      assert(len(matched) == len(predicted))
      assert(matched.count(True) <= len(reference))
      if matched.count(True) == 0:
        raise Sorry('''
          Invalid input for reference reflections.
          Zero reference spots were matched to predictions
        ''')
      elif matched.count(True) != len(reference):
        from logging import info
        info('')
        info('*' * 80)
        info('Warning: %d reference spots were not matched to predictions' % (
          len(reference) - matched.count(True)))
        info('*' * 80)
        info('')

      # Compute the profile model
      experiments = ProfileModelFactory.create(params, experiments, reference)

      # Compute the bounding box
      predicted.compute_bbox(experiments)

      # Create the integrator
      integrator = IntegratorFactory.create(params, experiments, predicted)

      # Integrate the reflections
      reflections = integrator.integrate()

      #print len(reflections)

      stats.integrated_intensity = flex.sum(reflections['intensity.sum.value'])

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
        response.append(
          '<fraction_indexed>%.2f</fraction_indexed>' %stats.fraction_indexed)
      if hasattr(stats, 'integrated_intensity'):
        response.append(
          '<integrated_intensity>%.0f</integrated_intensity>' %stats.integrated_intensity)

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

  # catch TERM signal to allow finalizers to run and reap daemonic children
  signal.signal(signal.SIGTERM, lambda *args: sys.exit(-signal.SIGTERM))

  for j in range(nproc - 1):
    proc = Process(target=serve, args=(httpd,))
    proc.daemon = True
    proc.start()
  serve(httpd)
  httpd.server_close()
  print time.asctime(), 'done'

if __name__ == '__main__':
  import sys
  import libtbx.load_env

  usage = '%s [options]' % libtbx.env.dispatcher_name

  from dials.util.options import OptionParser
  parser = OptionParser(
    usage=usage, phil=phil_scope, epilog=help_message)
  params, options = parser.parse_args(show_diff_phil=True)
  if params.nproc is libtbx.Auto:
    from libtbx.introspection import number_of_processors
    params.nproc = number_of_processors(return_value_if_unknown=-1)
  main(params.nproc, params.port)
