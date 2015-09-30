# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_tutorial_text

from __future__ import division
import os
import time
import shutil
import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx.test_utils import open_tmp_directory
from libtbx import easy_run

class Job(object):

  def __call__(self):
    print self.cmd
    result = self.timed_easy_run(self.cmd)
    print "{0:.3f}s".format(result.time)
    return result.stdout_lines

  @staticmethod
  def timed_easy_run(cmd):
    ts = time.time()
    result = easy_run.fully_buffered(command=cmd).raise_if_errors()
    te = time.time()
    result.time = te-ts
    return result

class dials_import(Job):

  def __init__(self):
    # find i04 bag training data, this may be part of dials_regression or xia2_regression
    if not libtbx.env.has_module("dials_regression") and not libtbx.env.has_module('xia2_regression'):
      raise RuntimeError("No dials_regression or xia2_regression module available!")

    data_dir = None

    # use the i04_weak_data for this test
    try:
      dials_regression = os.path.join(
        libtbx.env.dist_path('dials_regression'),
        "data", "i04-BAG-training")
      if os.path.isdir(dials_regression): data_dir = dials_regression
    except Exception:
      pass

    xia2_regression = os.path.join(
      abs(libtbx.env.build_path),
      'xia2_regression', 'test_data', 'i04_bag_training')
    if os.path.isdir(xia2_regression): data_dir = xia2_regression

    if data_dir is None:
      raise RuntimeError("Cannot find i04 data in either %s or %s" % (dials_regression, xia2_regression))
    self.cmd = "dials.import {0}".format(os.path.join(data_dir,"th_8_2_0*cbf"))

class dials_find_spots(Job):
  cmd = "dials.find_spots datablock.json nproc=4"

class dials_index(Job):
  cmd = "dials.index datablock.json strong.pickle"

class dials_refine_bravais_settings(Job):
  cmd = "dials.refine_bravais_settings experiments.json indexed.pickle"

class dials_refine(Job):
  cmd = "dials.refine bravais_setting_9.json indexed.pickle scan_varying=true"

class dials_integrate(Job):
  cmd = "dials.integrate refined_experiments.json refined.pickle background.algorithm=glm nproc=4"

class LogWriter(object):

  def __init__(self, log_dir):
    self.log_dir = log_dir

  def __call__(self, filename, text):
    with open(os.path.join(self.log_dir, filename), "w") as f:
      for line in text:
        f.write(line + "\n")

if (__name__ == "__main__") :

  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="generate_tutorial_text")
  os.chdir(tmp_dir)

  try:
    import_log = dials_import()()

    find_spots_log = dials_find_spots()()

    index_log = dials_index()()

    refine_bravais_settings_log = dials_refine_bravais_settings()()

    refine_log = dials_refine()()

    integrate_log = dials_integrate()()

    # if we got this far, assume it is okay to overwrite the logs
    dials_dir = libtbx.env.find_in_repositories("dials")
    log_dir = os.path.join(dials_dir, "doc", "sphinx", "documentation",
                           "tutorials", "logs")

    log_writer = LogWriter(log_dir)
    log_writer("dials.import.log", import_log)
    log_writer("dials.find_spots.log", find_spots_log)
    log_writer("dials.index.log", index_log)
    log_writer("dials.refine_bravais_settings.log", refine_bravais_settings_log)
    log_writer("dials.refine.log", refine_log)
    log_writer("dials.integrate.log", integrate_log)

    print "Updated log files written to {0}".format(log_dir)

  finally:
    os.chdir(cwd)
    # clean up tmp dir
    shutil.rmtree(tmp_dir)
