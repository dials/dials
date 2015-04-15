from __future__ import division
import os
import time
import libtbx.load_env
from libtbx import easy_run

have_dials_regression = libtbx.env.has_module("dials_regression")
if have_dials_regression:
  dials_regression = libtbx.env.find_in_repositories(
    relative_path="dials_regression",
    test=os.path.isdir)


def run():
  if not have_dials_regression:
    print "Skipping tst_find_spots_server_client: dials_regression not available."
    return

  server_command = "dials.find_spots_server"

  def start_server(server_command):
    result = easy_run.fully_buffered(command=server_command)
    result.show_stdout()
    result.show_stderr()

  import multiprocessing
  p = multiprocessing.Process(target=start_server, args=(server_command,))
  p.start()
  time.sleep(1) # need to give server chance to start

  try:
    exercise_client()

  finally:
    client_stop_command = "dials.find_spots_client stop"
    result = easy_run.fully_buffered(command=client_stop_command).raise_if_errors()
    #result.show_stdout()
    p.terminate()

def exercise_client():
  import glob
  data_dir = os.path.join(dials_regression, "centroid_test_data")
  filenames = glob.glob(os.path.join(data_dir, "*.cbf"))
  assert len(filenames) > 0
  client_command = " ".join(
    ["dials.find_spots_client",
     "min_spot_size=3",
     filenames[0]]
  )

  result = easy_run.fully_buffered(command=client_command).raise_if_errors()
  out = "<document>%s</document>" %"\n".join(result.stdout_lines)
  #result.show_stdout()

  from xml.dom import minidom
  xmldoc = minidom.parseString(out)
  assert len(xmldoc.getElementsByTagName('image')) == 1
  assert len(xmldoc.getElementsByTagName('spot_count')) == 1
  assert len(xmldoc.getElementsByTagName('spot_count_no_ice')) == 1
  assert len(xmldoc.getElementsByTagName('d_min')) == 1
  assert len(xmldoc.getElementsByTagName('total_intensity')) == 1

  client_command = " ".join([client_command] + filenames[1:])
  result = easy_run.fully_buffered(command=client_command).raise_if_errors()
  out = "<document>%s</document>" %"\n".join(result.stdout_lines)

  from xml.dom import minidom
  xmldoc = minidom.parseString(out)
  images = xmldoc.getElementsByTagName('image')
  assert len(images) == 9
  spot_counts = sorted([int(node.childNodes[0].data)
                 for node in xmldoc.getElementsByTagName('spot_count')])
  assert spot_counts == sorted([203, 196, 205, 209, 195, 205, 203, 207, 189]), spot_counts
  spot_counts_no_ice = sorted([
    int(node.childNodes[0].data)
    for node in xmldoc.getElementsByTagName('spot_count_no_ice')])
  assert spot_counts_no_ice \
         == sorted([150, 142, 151, 161, 151, 167, 164, 161, 146]), spot_counts_no_ice
  d_min = sorted([float(node.childNodes[0].data)
                  for node in xmldoc.getElementsByTagName('d_min')])
  assert d_min == sorted([1.47, 1.55, 1.59, 1.61, 1.61, 1.61, 1.61, 1.62, 1.64]), d_min


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    from libtbx.utils import show_times_at_exit
    show_times_at_exit()
    run()
