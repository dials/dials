# LIBTBX_SET_DISPATCHER_NAME dev.dials.semisynthetic_variance_analysis

from __future__ import division

def weighted_mean_variance(values, variances):
  import math
  weights = [1.0 / v for v in variances]
  sum_weights = sum(weights)
  weighted_mean = sum([v * w for v, w in zip(values, weights)]) / sum_weights
  weighted_variance = len(values) / sum_weights
  return weighted_mean, weighted_variance

def npp(values, input_mean_variance):
  import math
  from scitbx.math import distributions
  from scitbx.array_family import flex
  distribution = distributions.normal_distribution()
  values = flex.sorted(values)
  mean, variance = input_mean_variance
  scaled = (values - mean) / math.sqrt(variance)
  expected = distribution.quantiles(values.size())

  return expected, scaled

def semisynthetic_variance_analysis(semisynthetic_integrated_data_files,
                                    value_column):
  import cPickle as pickle
  import math
  from dials.array_family import flex
  from dials.util.add_hash import add_hash

  integrated_data_sets = [pickle.load(open(data_file, 'rb')) for
                          data_file in semisynthetic_integrated_data_files]

  # check column, find variance
  integrated_data = integrated_data_sets[0]
  assert value_column in integrated_data
  variance_column = None
  for column in ['%s.variance' % value_column,
                 value_column.replace('value', 'variance')]:
    if column in integrated_data:
      variance_column = column
      break
  assert(variance_column)
  data = integrated_data[value_column]
  if hasattr(data, 'parts'):
    multicolumn = len(data.parts())
  else:
    multicolumn = 0

  # first prepare the data files i.e. remove partials, keep only integrated
  # reflections, add the hash column, add weight column

  hash_set = None

  hashed_data_sets = []

  for integrated_data in integrated_data_sets:
    size0 = integrated_data.size()
    if 'intensity' in value_column:
      sel = integrated_data.get_flags(integrated_data.flags.integrated)
      integrated_data = integrated_data.select(sel)
      sel = integrated_data['partiality'] > 0.99
      integrated_data = integrated_data.select(sel)
    elif 'xyzobs' in value_column:
      sel = integrated_data.get_flags(integrated_data.flags.indexed)
      integrated_data = integrated_data.select(sel)
    integrated_data = add_hash(integrated_data)
    hashed_data_sets.append(integrated_data)
    if hash_set is None:
      hash_set = set(integrated_data['hash'])
    else:
      hash_set = hash_set.intersection(set(integrated_data['hash']))
    size1 = integrated_data.size()

  duplicate = []
  for h in hash_set:
    # check for duplicates i.e. reflection at 0, 2pi
    for i in hashed_data_sets:
      sel = i['hash'] == h
      isel = sel.iselection()
      if len(isel) > 1:
        duplicate.append(h)

  for d in duplicate:
    hash_set.discard(d)

  # now analyse those reflections found to be in all data sets (here looking
  # at the profile fitted intensity and variance thereof)

  for h in hash_set:
    if not multicolumn:
      values = flex.double()
      variances = flex.double()
      for i in hashed_data_sets:
        sel = i['hash'] == h
        isel = sel.iselection()
        assert(len(isel) == 1)
        values.append(i[isel[0]][value_column])
        variances.append(i[isel[0]][variance_column])
      weighted_mean, weighted_variance = weighted_mean_variance(values,
                                                                variances)
      expected, scaled = npp(values, (weighted_mean, weighted_variance))
      fit = flex.linear_regression(expected, scaled)
      # since I have everything needed to compute chi-square here...
      n = len(values)
      chi2 = sum([((v - weighted_mean) ** 2) / weighted_variance for v in values]) / n
      print '%.3f %.3f %.3f' % (weighted_mean / math.sqrt(weighted_variance), fit.slope(), chi2)
    else:
      values = { }
      variances = { }
      for m in range(multicolumn):
        values[m] = flex.double()
        variances[m] = flex.double()
      for i in hashed_data_sets:
        sel = i['hash'] == h
        isel = sel.iselection()
        assert(len(isel) == 1)
        data = i[isel[0]][value_column]
        variance = i[isel[0]][variance_column]
        for m in range(multicolumn):
          values[m].append(data[m])
          variances[m].append(variance[m])
      result = ''
      for m in range(multicolumn):
        weighted_mean, weighted_variance = weighted_mean_variance(values[m],
                                                                  variances[m])
        expected, scaled = npp(values[m], (weighted_mean, weighted_variance))
        fit = flex.linear_regression(expected, scaled)
        # since I have everything needed to compute chi-square here...
        n = len(values[m])
        chi2 = sum([((v - weighted_mean) ** 2) / weighted_variance for v in values[m]]) / n
        result += '%f %.3f %.3f ' % (math.sqrt(weighted_variance), fit.slope(), chi2)
      print result

if __name__ == '__main__':
  import sys
  value_column = sys.argv[1]
  semisynthetic_variance_analysis(sys.argv[2:], value_column)
