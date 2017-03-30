
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


def cluster_map(
    func,
    iterable,
    callback=None,
    nslots=1):
  '''
  A function to map stuff on cluster using drmaa

  :param func: The function to call
  :param iterable: The iterable to pass to each function call
  :param callback: A callback function when each job completes
  :param nslots: The number of processes to request per cluster node

  '''
  import multiprocessing
  import cPickle as pickle
  from os.path import join
  import os
  import sys
  import tempfile
  import drmaa

  # Set the working directory and make sure it exists
  # This will be where all the input/output files associated with the cluster
  # submission will go. For each job there will be a file:
  #  - ${INDEX}.input     The pickled input to the job
  #  - ${INDEX}.output    The pickled output from the job
  #  - ${INDEX}.stdout    The stdout from the job
  #  - ${INDEX}.stderr    The stderr from the job
  cwd = tempfile.mkdtemp(prefix="dials_cluster_map_", dir=os.getcwd())

  # Start outputting the input files in a separate process
  process = multiprocessing.Process(
    target = InputWriter(
      cwd,
      func,
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
    jt.outputPath       = ":" + join(cwd, "%s.stdout" % drmaa.JobTemplate.PARAMETRIC_INDEX)
    jt.errorPath        = ":" + join(cwd, "%s.stderr" % drmaa.JobTemplate.PARAMETRIC_INDEX)

    # FIXME Currently no portable way of specifying this
    # In order to select a cluster node with N cores
    # we have to use the native specification. This will work
    # on SGE but may not work for other queuing systems
    # This will set the NSLOTS environment variable
    jt.nativeSpecification = '-pe smp %d' % nslots

    N = len(iterable)
    try:

      # Submit the array job
      joblist = s.runBulkJobs(jt, 1, N, 1)

      # For each item, load the result and process
      result = []
      for i, jobid in enumerate(joblist, start=1):
        s.wait(jobid, drmaa.Session.TIMEOUT_WAIT_FOREVER)
        with open(join(cwd, "%d.stdout" % i), "r") as infile:
          sys.stdout.write("".join(infile.readlines()))
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
  from dials.util.mp import MultiNodeClusterFunction



  print cluster_map(
    func,
    list(range(100)),
    nslots=4
   # callback=callback,
    )
