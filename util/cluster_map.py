
from __future__ import absolute_import, division

class InputWriter(object):
  '''
  A class to write the input files

  '''

  def __init__(self, directory, function, iterable):
    '''
    Save the function and iterable

    '''
    self.directory = directory
    self.function = function
    self.iterable = iterable

  def __call__(self):
    '''
    Call this to write input files

    '''
    from os.path import join
    import cPickle as pickle
    for i, item in enumerate(self.iterable, start=1):
      with open(join(self.directory, "%d.input" % i), "wb") as outfile:
        pickle.dump(
          (self.function, item),
          outfile,
          protocol=pickle.HIGHEST_PROTOCOL)


def cluster_map(function, iterable, callback=None):
  '''
  A function to map stuff on cluster using drmaa

  '''
  import multiprocessing
  import cPickle as pickle
  from os.path import join
  import os
  import drmaa

  # Set the working directory and make sure it exists
  cwd = join(os.getcwd(), "cluster_io_data")
  try:
    os.mkdir(cwd)
  except OSError as exc:
    import errno
    if exc.errno != errno.EEXIST:
      raise exc

  # Start outputting the input files in a separate process
  process = multiprocessing.Process(
    target = InputWriter(
      cwd,
      function,
      iterable))
  process.start()

  # Start the drmaa session
  with drmaa.Session() as s:

    # Create the job template
    jt = s.createJobTemplate()
    jt.remoteCommand    = "cluster.dials.exec"
    jt.args             = [cwd]
    jt.jobName          = "dials"
    jt.joinFiles        = True
    jt.jobEnvironment   = os.environ
    jt.workingDirectory = cwd
    #jt.nativeSpecification = '-pe smp %d' % 10

    N = len(iterable)
    try:
      joblist = s.runBulkJobs(jt, 1, N, 1)

      # For each item, load the result and process
      result = []
      for i, jobid in enumerate(joblist, start=1):
        s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        with open(join(cwd, "%d.output" % i), "rb") as infile:
          r = pickle.load(infile)
        if isinstance(r, Exception):
          raise r
        if callback is not None:
          callback(r)
        result.append(r)

      # Make sure the process has finished
      process.join()

      # Delete job template
      s.deleteJobTemplate(jt)

    except KeyboardInterrupt:

      # Delete the jobs
      s.control(
        drmaa.Session.JOB_IDS_SESSION_ALL,
        drmaa.JobControlAction.TERMINATE)

      # Delete job template
      s.deleteJobTemplate(jt)

      # Re-raise exception
      raise

  # Return the result
  return result

if __name__ == '__main__':

  def callback(x):
    print x

  from dials.util.cluster_func_test import func
  print cluster_map(func, list(range(10000)), callback=callback)
