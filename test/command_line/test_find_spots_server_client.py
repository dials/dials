from __future__ import absolute_import, division, print_function

import glob
import multiprocessing
import os
import procrunner
import pytest
import socket
import time
import timeit
from xml.dom import minidom

from libtbx import easy_run

def start_server(server_command):
  procrunner.run(server_command)

def test_find_spots_server_client(dials_regression, tmpdir):
  tmpdir.chdir()

  s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
  s.bind(("",0))
  port = s.getsockname()[1]
  server_command = ["dials.find_spots_server", "port=%i" % port, "nproc=3"]
  print(server_command)

  p = multiprocessing.Process(target=start_server, args=(server_command,))
  p.daemon = True
  s.close()
  p.start()
  wait_for_server(port) # need to give server chance to start

  data_dir = os.path.join(dials_regression, "centroid_test_data")
  filenames = sorted(glob.glob(os.path.join(data_dir, "*.cbf")))

  try:
    exercise_client(port=port, filenames=filenames)

  finally:
    client_stop_command = "dials.find_spots_client port=%i stop" %port
    result = easy_run.fully_buffered(command=client_stop_command).raise_if_errors()
    #result.show_stdout()
    p.terminate()

def wait_for_server(port, max_wait=20):
  print("Waiting up to %d seconds for server to start" % max_wait)
  server_ok = False
  start_time = timeit.default_timer()
  max_time = start_time + max_wait
  while (timeit.default_timer() < max_time) and not server_ok:
    try:
      s = socket.socket(socket.AF_INET, socket.SOCK_STREAM)
      s.connect(('127.0.0.1', port))
      s.close()
      server_ok = True
    except socket.error as e:
      if (e.errno != 111) and (e.errno != 61):
        raise
      # ignore connection failures (111 connection refused on linux; 61 connection refused on mac)
      time.sleep(0.1)
  if not server_ok:
    raise Exception("Server failed to start after %d seconds" % max_wait)
  print("dials.find_spots_server up after %f seconds" % (timeit.default_timer() - start_time))

def exercise_client(port, filenames):
  assert len(filenames) > 0
  client_command = " ".join(
    ["dials.find_spots_client",
     "port=%i" %port,
     "min_spot_size=3",
     "nproc=1",
     filenames[0]]
  )

  index_client_command = " ".join(
    [client_command,
     "index=True",
     "indexing.method=fft1d",
     "max_refine=10",])
  print(index_client_command)
  result = easy_run.fully_buffered(command=index_client_command).raise_if_errors()
  out = "<document>%s</document>" %"\n".join(result.stdout_lines)
  result.show_stdout()
  #result.show_stderr()

  xmldoc = minidom.parseString(out)
  assert len(xmldoc.getElementsByTagName('image')) == 1
  assert len(xmldoc.getElementsByTagName('spot_count')) == 1
  assert len(xmldoc.getElementsByTagName('spot_count_no_ice')) == 1
  assert len(xmldoc.getElementsByTagName('d_min')) == 1
  assert len(xmldoc.getElementsByTagName('total_intensity')) == 1
  assert len(xmldoc.getElementsByTagName('unit_cell')) == 1
  assert len(xmldoc.getElementsByTagName('n_indexed')) == 1
  assert len(xmldoc.getElementsByTagName('fraction_indexed')) == 1

  unit_cell = [
    float(f) for f in
    xmldoc.getElementsByTagName('unit_cell')[0].childNodes[0].data.split()]

  assert unit_cell == pytest.approx(
    [39.90, 42.67, 42.37, 89.89, 90.10, 90.13], abs=1e-1)

  client_command = " ".join([client_command] + filenames[1:])
  result = easy_run.fully_buffered(command=client_command).raise_if_errors()
  out = "<document>%s</document>" %"\n".join(result.stdout_lines)

  xmldoc = minidom.parseString(out)
  images = xmldoc.getElementsByTagName('image')
  assert len(images) == 9
  spot_counts = sorted([int(node.childNodes[0].data)
                 for node in xmldoc.getElementsByTagName('spot_count')])
  assert spot_counts == sorted([203, 196, 205, 209, 195, 205, 203, 207, 189])
  spot_counts_no_ice = sorted([
    int(node.childNodes[0].data)
    for node in xmldoc.getElementsByTagName('spot_count_no_ice')])
  assert spot_counts_no_ice \
         == sorted([169, 171, 175, 176, 177, 184, 193, 195, 196])
  d_min = sorted([float(node.childNodes[0].data)
                  for node in xmldoc.getElementsByTagName('d_min')])
  assert d_min == sorted([1.45, 1.47, 1.55, 1.55, 1.56, 1.59, 1.61, 1.61, 1.64])
