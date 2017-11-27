# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_tutorial_text

from __future__ import absolute_import, division
import sys
import os
import glob
import shutil
import shlex
import argparse

import libtbx.load_env # required for libtbx.env.find_in_repositories
from libtbx.test_utils import open_tmp_directory
from dials.util.procrunner import run_process

class Job(object):
  """Represents a step command to execute"""

  def __call__(self):
    """Run the command this job represents.

    Standard output is saved onto self.result, and then optionally filtered
    by the mangle_result function.

    :returns: Dictionary with keys 'cmd' and 'result'
    """
    print self.cmd
    self.result = Job.run_process(self.cmd)
    self.mangle_result()
    return {'cmd':self.cmd, 'result':self.result['stdout']}

  def mangle_result(self):
    ''' function that can be overridden to change the return values after execution '''
    pass

  @staticmethod
  def run_process(command):
    """Runs a command, prints running info and results the result, if success"""
    result = run_process(shlex.split(command))
    print "running command took {0:.2f} seconds\n".format(result['runtime'])
    assert result['exitcode'] == 0, "Command execution failed"
    return result

class JobWriter(object):
  '''Tool to save job command and result to files in a fixed destination'''

  def __init__(self, directory):
    self.directory = directory

  def __call__(self, cmd_filename, result_filename, job):
    """
    Save command and output to named files.

    :param cmd_filename:    Where to save a copy of the run command
    :param result_filenam:  Where to save the result output
    :param job:             The result dictionary from calling a Job
    """
    with open(os.path.join(self.directory, cmd_filename), "w") as f:
      f.write(job['cmd'])
    with open(os.path.join(self.directory, result_filename), "w") as f:
      f.write(job['result'])


class Processing_Tutorial(object):
  """Command steps for generating the logs for the processing tutorial"""

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

  class dials_reindex(Job):
    cmd = "dials.reindex indexed.pickle change_of_basis_op=a,b,c"

  class dials_refine(Job):
    cmd = "dials.refine bravais_setting_9.json reindexed_reflections.pickle"

  class dials_sv_refine(Job):
    cmd = "dials.refine refined_experiments.json refined.pickle scan_varying=true"

  class dials_integrate(Job):
    cmd = "dials.integrate refined_experiments.json refined.pickle nproc=4"

  class dials_report(Job):
    cmd = "dials.report integrated_experiments.json integrated.pickle"

    def mangle_result(self):
      self.result['stdout'] = open('dials-report.html').read()

  class dials_export(Job):
    cmd = "dials.export integrated.pickle integrated_experiments.json mtz.hklout=integrated.mtz"

def generate_processing_detail_text():
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="generate_tutorial_text")
  os.chdir(tmp_dir)

  try:
    import_job = dials_import()()

    find_spots_job = dials_find_spots()()

    index_job = dials_index()()

    refine_bravais_settings_job = dials_refine_bravais_settings()()

    reindex_job = dials_reindex()()

    refine_job = dials_refine()()

    sv_refine_job = dials_sv_refine()()

    integrate_job = dials_integrate()()

    report_html_job = dials_report()()

    export_job = dials_export()()

    # if we got this far, assume it is okay to overwrite the logs
    dials_dir = libtbx.env.find_in_repositories("dials")
    result_dir = os.path.join(dials_dir, "doc", "sphinx", "documentation",
                           "tutorials", "logs")

    job_writer = JobWriter(result_dir)
    job_writer("dials.import.cmd", "dials.import.log", import_job)
    job_writer("dials.find_spots.cmd", "dials.find_spots.log", find_spots_job)
    job_writer("dials.index.cmd", "dials.index.log", index_job)
    job_writer("dials.refine_bravais_settings.cmd",
               "dials.refine_bravais_settings.log", refine_bravais_settings_job)
    job_writer("dials.reindex.cmd", "dials.reindex.log", reindex_job)
    job_writer("dials.refine.cmd","dials.refine.log", refine_job)
    job_writer("dials.sv_refine.cmd","dials.sv_refine.log", sv_refine_job)
    job_writer("dials.integrate.cmd", "dials.integrate.log", integrate_job)
    job_writer("dials.report.cmd", "dials-report.html", report_html_job)
    job_writer("dials.export.cmd", "dials.export.log", export_job)

    print "Updated result files written to {0}".format(result_dir)

  finally:
    os.chdir(cwd)
    # clean up tmp dir
    shutil.rmtree(tmp_dir)

def generate_processing_detail_text_ccp4():
  """Generate the text for CCP4-versions of detail processing tutorial"""

  # Move to a temporary directory for processing
  cwd = os.path.abspath(os.curdir)
  tmp_dir = open_tmp_directory(suffix="generate_tutorial_text_ccp4")
  os.chdir(tmp_dir)

  # Find/validate the data input - until we've decided to integrate this
  # into the main release, have a DLS default or otherwise let it be
  # specified via a CCP4_TUTORIAL_DATA environment variable.
  DATA_PATH = "/dls/i03/data/2017/mx19576-1/tutorial_data/summed/C2sum_1*.cbf.gz"
  DATA_PATH = os.environ.get("CCP4_TUTORIAL_DATA", DATA_PATH)
  if not any(glob.glob(DATA_PATH)):
    print("Error: Could not find CCP4-2017 data; skipping text generation")
    return

  # Work out where we are writing the output files to; in-source
  dials_dir = libtbx.env.find_in_repositories("dials")
  OUTPUT_DIR = os.path.join(dials_dir, "doc", "sphinx", "documentation",
                            "tutorials", "logs_ccp4")
  # Ensure this output path exists
  if not os.path.isdir(OUTPUT_DIR):
    os.makedirs(OUTPUT_DIR)

  # Make an ordered list of named steps and associated commands
  commands = [
    ("dials.import",                  "dials.import {}".format(DATA_PATH)),
    ("dials.find_spots",              "dials.find_spots datablock.json nproc=4"),
    ("dials.index",                   "dials.index datablock.json strong.pickle"),
    ("dials.refine_bravais_settings", "dials.refine_bravais_settings experiments.json indexed.pickle"),
    ("dials.reindex",                 "dials.reindex indexed.pickle change_of_basis_op=a+b,-a+b,c"),
    ("dials.refine",                  "dials.refine bravais_setting_2.json reindexed_reflections.pickle"),
    ("dials.sv_refine",               "dials.refine refined_experiments.json refined.pickle scan_varying=true"),
    ("dials.integrate",               "dials.integrate refined_experiments.json refined.pickle nproc=4"),
    ("dials.report",                  "dials.report integrated_experiments.json integrated.pickle"),
    ("dials.export",                  "dials.export integrated.pickle integrated_experiments.json mtz.hklout=integrated.mtz"),
  ]

  job_writer = JobWriter(OUTPUT_DIR)

  # Protect against errors so we can tidy up afterwards
  try:
    # Run each step, and write an output log
    for name, command in commands:
      result = Job.run_process(command)["stdout"]
      # Write a copy of the command, and the output log
      job_writer("{}.cmd".format(name), "{}.log".format(name),{'cmd':command,'result': result})

    # Report step is special; we want the dials-report.html file instead
    shutil.copy("dials-report.html", OUTPUT_DIR)

    print "Updated result files written to {0}".format(OUTPUT_DIR)

  finally:
    # Remove our intermediatary files
    os.chdir(cwd)
    shutil.rmtree(tmp_dir)

if __name__ == "__main__":
  # As a quick development hack, add option for only the newer process
  if not "--ccp4" in sys.argv:
    generate_processing_detail_text()
  generate_processing_detail_text_ccp4()
