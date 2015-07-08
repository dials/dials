from logging import info, debug
from libtbx.phil import parse
from dials.array_family import flex

class CentroidOutlier(object):
  """Base class for centroid outlier detection algorithms"""

  def __init__(self, cols=["x_resid", "y_resid", "phi_resid"],
               min_num_obs=20,
               separate_experiments=True,
               separate_panels=True):

    # column names of the data in which to look for outliers
    self._cols = cols

    # minimum number of observations per panel below which all reflections will
    # be marked as potential outliers
    self._min_num_obs = min_num_obs

    # whether to do outlier rejection within each experiment or panel, or across
    # all experiments and panels
    self._separate_experiments = separate_experiments
    self._separate_panels = separate_panels

    # the number of rejections
    self.nreject = 0

    return

  def _detect_outliers(cols):
    """Perform outlier detection using the input cols and return a flex.bool
    indicating which rows in the cols are considered outlying. cols should be
    a list of flex.doubles of equal lengths"""

    # to be implemented by derived classes
    raise NotImplementedError()

  def __call__(self, reflections):
    """Identify outliers in the input and set the centroid_outlier flag.
    Return True if any outliers were detected, otherwise False"""

    # check the columns are present
    for col in self._cols: assert reflections.has_key(col)

    sel = reflections.get_flags(reflections.flags.used_in_refinement)
    all_data = reflections.select(sel)
    all_data_indices = sel.iselection()

    jobs = []
    if self._separate_experiments:
      # split the data set by experiment id
      for iexp in xrange(flex.max(all_data['id'] + 1)):
        sel = all_data['id'] == iexp
        job = {'id':iexp, 'panel':'all', 'data':all_data.select(sel),
               'indices':all_data_indices.select(sel)}
        jobs.append(job)
    else:
      # keep the whole dataset across all experiment ids
      job = {'id':'all', 'panel':'all', 'data':all_data,
             'indices':all_data_indices}
      jobs.append(job)

    jobs2 = []
    if self._separate_panels:
      # split further by panel id
      for job in jobs:
        data = job['data']
        iexp = job['id']
        indices = job['indices']
        for ipanel in xrange(flex.max(data['panel']) + 1):
          sel = data['panel'] == ipanel
          job2 = {'id':iexp, 'panel':ipanel, 'data':data.select(sel),
                  'indices':indices.select(sel)}
          jobs2.append(job2)
    else:
      # keep the splits as they are
      jobs2 = jobs

    # now loop over the lowest level of splits
    for job in jobs2:

      data = job['data']
      indices = job['indices']
      iexp = job['id']
      ipanel = job['panel']

      if len(indices) >= self._min_num_obs:

        # get the subset of data as a list of columns
        cols = [data[col] for col in self._cols]

        # determine the position of outliers on this sub-dataset
        outliers = self._detect_outliers(cols)

        # get positions of outliers from the original matches
        ioutliers = indices.select(outliers)

      else:
        msg = "For experiment: {0} and panel: {1}, ".format(iexp, ipanel)
        msg += "only {0} reflections are present. ".format(len(indices))
        msg += "All of these flagged as possible outliers."
        debug(msg)
        ioutliers = indices

      # set those reflections as outliers in the original reflection table
      reflections.set_flags(ioutliers,
        reflections.flags.centroid_outlier)

      self.nreject += len(ioutliers)

    if self.nreject == 0: return False

    info("{0} reflections have been flagged as outliers".format(self.nreject))

    return True

# The phil scope for outlier rejection
phil_str = '''
outlier
  .help = "Outlier rejection after initial reflection prediction."
{
  algorithm = *null mcd tukey
    .help = "Outlier rejection algorithm"
    .type = choice

  minimum_number_of_reflections = 20
    .help = "The minimum number of input observations per outlier rejection"
            "job below which all reflections in the job will be rejected as"
            "potential outliers"
    .type = int(value_min=0)

  separate_experiments = True
    .help = ""
    .type = bool

  separate_panels = True
    .help = ""
    .type = bool

  tukey
    .help = "Options for the tukey outlier rejector"
    .expert_level = 1
  {
    iqr_multiplier = 1.5
      .help = "The IQR multiplier used to detect outliers. A value of 1.5"
              "gives Tukey's rule for outlier detection"
      .type = float(value_min = 0.)
  }

  mcd
    .help = "Options for the mcd outlier rejector"
    .expert_level = 1
  {
     alpha = 0.5
       .help = ""
       .type = float(value_min = 0., value_max = 1.0)

     max_n_groups=5
       .help = ""
       .type = int(value_min = 1)

     min_group_size=300
       .help = ""
       .type = int(value_min = 100)

     n_trials=500
       .help = ""
       .type = int(value_min = 1)

     k1=2
       .help = ""
       .type = int(value_min = 1)

     k2=2
       .help = ""
       .type = int(value_min = 1)

     k3=100
       .help = ""
       .type = int(value_min = 1)

     threshold_probability=0.975
       .help = ""
       .type = float(value_min = 0., value_max = 1.0)
  }
}

'''

phil_scope = parse(phil_str)

class CentroidOutlierFactory(object):

  @classmethod
  def from_parameters_and_colnames(cls, params, colnames):

    # id the relevant scope for the requested method
    method = params.outlier.algorithm
    if method == "null": return None
    elif method == "tukey":
      from tukey import Tukey as outlier_detector
      algo_params = params.outlier.tukey
    elif method == "mcd":
      from mcd import MCD as outlier_detector
      algo_params = params.outlier.mcd
    else:
      raise RuntimeError("outlier.algorithm not recognised")

    # construct kwargs from the algo_params scope
    kwargs = dict((k, v) for k, v in algo_params.__dict__.items() \
      if not k.startswith('_'))

    return outlier_detector(
      cols=colnames,
      min_num_obs=params.outlier.minimum_number_of_reflections,
      separate_experiments=params.outlier.separate_experiments,
      separate_panels=params.outlier.separate_panels,
      **kwargs)

if __name__ == "__main__":

  params=phil_scope.extract()
  params.outlier.algorithm="tukey"
  print CentroidOutlierFactory.from_parameters_and_colnames(params, [1, 2, 3])

