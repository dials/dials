#
# report.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division
from dials.array_family import flex
from dials.array_family.flex import Binner


def flex_ios(val, var):
  '''
  Compute I/sigma or return zero for each element.

  '''
  assert(len(val) == len(var))
  result = flex.double(len(val),0)
  indices = flex.size_t(range(len(val))).select(var > 0)
  val = val.select(indices)
  var = var.select(indices)
  assert(var.all_gt(0))
  result.set_selected(indices, val / flex.sqrt(var))
  return result


def generate_integration_report(experiment, reflections, n_resolution_bins=20):
  '''
  Generate the integration report

  '''
  from collections import OrderedDict
  from dials.algorithms.statistics import pearson_correlation_coefficient
  from dials.algorithms.statistics import spearman_correlation_coefficient
  from cctbx import miller, crystal

  def overall_report(data):

    # Start by adding some overall numbers
    report = OrderedDict()
    report['n']           = len(reflections)
    report['n_full']      = data['full'].count(True)
    report['n_partial']   = data['full'].count(False)
    report['n_overload']  = data['over'].count(True)
    report['n_ice']       = data['ice'].count(True)
    report['n_summed']    = data['sum'].count(True)
    report['n_fitted']    = data['prf'].count(True)
    report['n_integated'] = data['int'].count(True)

    # Compute mean background
    try:
      report['mean_background'] = flex.mean(
        data['background.mean'].select(data['int']))
    except Exception:
      report['mean_background'] = 0.0

    # Compute mean I/Sigma summation
    try:
      report['ios_sum'] = flex.mean(
        data['intensity.sum.ios'].select(data['sum']))
    except Exception:
      report['ios_sum'] = 0.0

    # Compute mean I/Sigma profile fitting
    try:
      report['ios_prf'] = flex.mean(
        data['intensity.prf.ios'].select(data['prf']))
    except Exception:
      report['ios_prf'] = 0.0

    # Compute the mean profile correlation
    try:
      report['cc_prf'] = flex.mean(
        data['profile.correlation'].select(data['prf']))
    except Exception:
      report['cc_prf'] = 0.0

    # Compute the correlations between summation and profile fitting
    try:
      mask = data['sum'] & data['prf']
      Isum = data['intensity.sum.value'].select(mask)
      Iprf = data['intensity.prf.value'].select(mask)
      report['cc_pearson_sum_prf'] = pearson_correlation_coefficient(Isum, Iprf)
      report['cc_spearman_sum_prf'] = spearman_correlation_coefficient(Isum, Iprf)
    except Exception:
      report['cc_pearson_sum_prf'] = 0.0
      report['cc_spearman_sum_prf'] = 0.0

    # Return the overall report
    return report

  def binned_report(binner, index, data):

    # Create the indexers
    indexer_all = binner.indexer(index)
    indexer_sum = binner.indexer(index.select(data['sum']))
    indexer_prf = binner.indexer(index.select(data['prf']))
    indexer_int = binner.indexer(index.select(data['int']))

    # Add some stats by resolution
    report = OrderedDict()
    report['bins']         = list(binner.bins())
    report['n_full']       = list(indexer_all.sum(data['full']))
    report['n_partial']    = list(indexer_all.sum(~data['full']))
    report['n_overload']   = list(indexer_all.sum(data['over']))
    report['n_ice']        = list(indexer_all.sum(data['ice']))
    report['n_summed']     = list(indexer_all.sum(data['sum']))
    report['n_fitted']     = list(indexer_all.sum(data['prf']))
    report['n_integrated'] = list(indexer_all.sum(data['int']))

    # Compute mean background
    try:
      report['mean_background'] = list(indexer_int.mean(
        data['background.mean'].select(data['int'])))
    except Exception:
      report['mean_background'] = [0.0] * len(binner)

    # Compute mean I/Sigma summation
    try:
      report['ios_sum'] = list(indexer_sum.mean(
        data['intensity.sum.ios'].select(data['sum'])))
    except Exception:
      report['ios_sum'] = [0.0] * len(binner)

    # Compute mean I/Sigma profile fitting
    try:
      report['ios_prf'] = list(indexer_prf.mean(
        data['intensity.prf.ios'].select(data['prf'])))
    except Exception:
      report['ios_prf'] = [0.0] * len(binner)

    # Compute the mean profile correlation
    try:
      report['cc_prf'] = list(indexer_prf.mean(
        data['profile.correlation'].select(data['prf'])))
    except Exception:
      report['cc_prf'] = [0.0] * len(binner)

    # Return the binned report
    return report

  def resolution_bins(experiment, hkl, nbins):

    # Create the crystal symmetry object
    cs = crystal.symmetry(
      space_group=experiment.crystal.get_space_group(),
      unit_cell=experiment.crystal.get_unit_cell())

    # Create the resolution binner object
    ms = miller.set(cs, hkl)
    ms.setup_binner(n_bins=nbins)
    binner = ms.binner()
    brange = list(binner.range_used())
    bins = [binner.bin_d_range(brange[0])[0]]
    for i in brange:
      bins.append(binner.bin_d_range(i)[1])
    return flex.double(reversed(bins))

  def select(data, indices):

    # Select rows from columns
    result = {}
    for key, value in data.iteritems():
      result[key] = value.select(indices)
    return result

  # Check the required columns are there
  assert("miller_index" in reflections)
  assert("d" in reflections)
  assert("flags" in reflections)
  assert("bbox" in reflections)
  assert("xyzcal.px" in reflections)
  assert("partiality" in reflections)
  assert("intensity.sum.value" in reflections)
  assert("intensity.sum.variance" in reflections)

  # Get the flag enumeration
  flags = flex.reflection_table.flags

  # Get some keys from the data
  data = {}
  for key in ['miller_index',
              'xyzcal.px',
              'd',
              'bbox',
              'background.mean',
              'partiality',
              'intensity.sum.value',
              'intensity.sum.variance',
              'intensity.prf.value',
              'intensity.prf.variance',
              'profile.correlation']:
    if key in reflections:
      data[key] = reflections[key]

  # Compute some flag stuff
  data["full"] = data['partiality'] > 0.997300203937
  data["over"] = reflections.get_flags(flags.overloaded)
  data["ice"] = reflections.get_flags(flags.in_powder_ring)
  data["sum"] = reflections.get_flags(flags.integrated_sum)
  data["prf"] = reflections.get_flags(flags.integrated_prf)
  data["int"] = reflections.get_flags(flags.integrated, all=False)

  # Try to calculate the i over sigma for summation
  data['intensity.sum.ios'] = flex_ios(
    data['intensity.sum.value'],
    data['intensity.sum.variance'])

  # Try to calculate the i over sigma for profile fitting
  try:
    data['intensity.prf.ios'] = flex_ios(
      data['intensity.prf.value'],
      data['intensity.prf.variance'])
  except Exception:
    pass

  # Create the resolution binner
  high_low_resolution_binner = Binner(
    resolution_bins(
      experiment,
      data['miller_index'],
      2))

  # Create the resolution binner
  resolution_binner = Binner(
    resolution_bins(
      experiment,
      data['miller_index'],
      n_resolution_bins))

  # Create the frame binner object
  try:
    array_range = experiment.imageset.get_array_range()
  except:
    array_range = (0, len(experiment.imageset))
  frame_binner = Binner(flex.int(range(
    array_range[0],
    array_range[1]+1)).as_double())

  # Create the overall report
  overall = overall_report(data)

  # Create high/low resolution reports
  hl_binner = high_low_resolution_binner.indexer(data['d'])
  high_summary = overall_report(select(data, hl_binner.indices(0)))
  low_summary = overall_report(select(data, hl_binner.indices(1)))

  # Create the overall report
  summary = OrderedDict([
    ('overall', overall),
    ('low', low_summary),
    ('high', high_summary)
  ])

  # Create a report binned by resolution
  resolution = binned_report(resolution_binner, data['d'], data)

  # Create the report binned by image
  image = binned_report(frame_binner, data['xyzcal.px'].parts()[2], data)

  # Return the report
  return OrderedDict([
    ("summary", summary),
    ("resolution", resolution),
    ("image", image)])


class IntegrationReport(object):
  '''
  A class to store the integration report

  '''

  def __init__(self, experiments, reflections):
    '''
    Create the integration report

    :param experiments: The experiment list
    :param reflections: The reflection table

    '''
    from collections import OrderedDict

    # Split the tables by experiment id
    tables = reflections.split_by_experiment_id()
    assert(len(tables) == len(experiments))

    # Initialise the dictionary
    self._report = []

    # Generate an integration report for each experiment
    for i, (expr, data) in enumerate(zip(experiments, tables)):
      self._report.append(generate_integration_report(expr, data))

  def as_dict(self):
    '''
    Return the report as a dictionary

    :return: The report dictionary

    '''
    return self._report

  def as_str(self, prefix=''):
    '''
    Return the report as a string

    :return: The report string

    '''
    from libtbx.table_utils import format as table

    # Create the image table
    rows = [["Id",
             "Image",
             "# full",
             "# part",
             "# over",
             "# ice",
             "# sum",
             "# prf",
             "<Ibg>",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)",
             "<CC prf>"]]
    for j, report in enumerate(self._report):
      report = report['image']
      for i in range(len(report['bins'])-1):
        rows.append([
          '%d'   % j,
          '%d'   % report['bins'][i],
          '%d'   % report['n_full'][i],
          '%d'   % report['n_partial'][i],
          '%d'   % report['n_overload'][i],
          '%d'   % report['n_ice'][i],
          '%d'   % report['n_summed'][i],
          '%d'   % report['n_fitted'][i],
          '%.2f' % report['mean_background'][i],
          '%.2f' % report['ios_sum'][i],
          '%.2f' % report['ios_prf'][i],
          '%.2f' % report['cc_prf'][i]])
    image_table = table(rows, has_header=True, justify='right', prefix=prefix)

    # Create the resolution table
    rows = [["Id",
             "d min",
             "# full",
             "# part",
             "# over",
             "# ice",
             "# sum",
             "# prf",
             "<Ibg>",
             "<I/sigI>\n (sum)",
             "<I/sigI>\n (prf)",
             "<CC prf>"]]
    for j, report in enumerate(self._report):
      report = report['resolution']
      for i in range(len(report['bins'])-1):
        rows.append([
          '%d'   % j,
          '%.2f' % report['bins'][i],
          '%d'   % report['n_full'][i],
          '%d'   % report['n_partial'][i],
          '%d'   % report['n_overload'][i],
          '%d'   % report['n_ice'][i],
          '%d'   % report['n_summed'][i],
          '%d'   % report['n_fitted'][i],
          '%.2f' % report['mean_background'][i],
          '%.2f' % report['ios_sum'][i],
          '%.2f' % report['ios_prf'][i],
          '%.2f' % report['cc_prf'][i]])
    resolution_table = table(rows, has_header=True, justify='right', prefix=prefix)

    # Create the overall table
    overall_tables = []
    for j, report in enumerate(self._report):
      report = report['summary']
      summary = report['overall']
      high = report['high']
      low = report['low']
      desc_fmt_key = [
        ("number fully recorded",                 '%d'  , "n_full"),
        ("number partially recorded",             '%d'  , "n_partial"),
        ("number with overloaded pixels",         '%d'  , "n_overload"),
        ("number in powder rings",                '%d'  , "n_ice"),
        ("number processed with summation",       '%d'  , "n_summed"),
        ("number processed with profile fitting", '%d'  , "n_fitted"),
        ("<ibg>",                                 '%.2f', "mean_background"),
        ("<i/sigi> (summation)",                  '%.2f', "ios_sum"),
        ("<i/sigi> (profile fitting)",            '%.2f', "ios_prf"),
        ("<cc prf>",                              '%.2f', "cc_prf"),
        ("cc_pearson sum/prf",                    '%.2f', "cc_pearson_sum_prf"),
        ("cc_spearman sum/prf",                   '%.2f', "cc_spearman_sum_prf")
      ]
      rows = [['item', 'overall', 'low', 'high']]
      for desc, fmt, key in desc_fmt_key:
        rows.append([desc, fmt % summary[key], fmt % low[key], fmt % high[key]])
      overall_tables.append((j, table(rows, has_header=True, justify='left', prefix=prefix)))

    # Create the text
    text = [
      prefix + 'Summary vs image number',
      image_table,
      '\n',
      prefix + 'Summary vs resolution',
      resolution_table,
      '\n']
    for i, table in overall_tables:
      text.append(prefix + 'Summary for experiment %d' % i)
      text.append(table)
      text.append('\n')

    # Return the text
    return '\n'.join(text)


class ProfileModelReport(object):
  '''
  A class to store the profile model report

  '''

  def __init__(self, experiments, profile_model, reflections):
    '''
    Create the integration report

    :param experiments: The experiment list
    :param profile_model: The profile model
    :param reflections: The reflection table

    '''
    from collections import OrderedDict

    # Init the report
    self._report = []

    # Get the modeller
    profiles = profile_model.profiles()

    # Create the summary for each profile model
    for i in range(len(profiles)):
      model = profiles[i]
      summary = OrderedDict()
      summary['valid'] = [model.valid(i) for i in range(len(model))]
      summary['coord'] = [model.coord(i) for i in range(len(model))]
      summary['n_reflections'] = [model.n_reflections(i) for i in range(len(model))]
      self._report.append(summary)

  def as_dict(self):
    '''
    Return the report as a dictionary

    :return: The report dictionary

    '''
    return self._report

  def as_str(self, prefix=''):
    '''
    Return the report as a string

    :return: The report string

    '''
    from libtbx.table_utils import format as table

    # Create the image table
    rows = [["Id",
             "Profile",
             "Created",
             "X (px)",
             "Y (px)",
             "Z (im)",
             "# reflections"]]
    for j, report in enumerate(self._report):
      for i in range(len(report['valid'])):
        rows.append([
          '%d'   % j,
          '%d'   % i,
          '%r'   % report['valid'][i],
          '%.2f' % report['coord'][i][0],
          '%.2f' % report['coord'][i][1],
          '%.2f' % report['coord'][i][2],
          '%d'   % report['n_reflections'][i]])
    summary_table = table(rows, has_header=True, justify='right', prefix=prefix)

    # Create the text
    text = [
      'Summary of profile model',
      summary_table
    ]

    # Return the text string
    return '\n'.join(text)
