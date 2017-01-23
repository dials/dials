# LIBTBX_SET_DISPATCHER_NAME dev.dials.saturation_analysis

from __future__ import division

def strip_not_integrated(integrated_data):
  sel = integrated_data.get_flags(integrated_data.flags.integrated)
  integrated_data = integrated_data.select(sel)
  sel = integrated_data['partiality'] > 0.99
  integrated_data = integrated_data.select(sel)
  return integrated_data

def saturation_analysis(data_files, value_column):
  import cPickle as pickle
  import math
  from dials.array_family import flex
  from dials.util.add_hash import add_hash, dehash
  from annlib_ext import AnnAdaptor as ann_adaptor

  reference = data_files[0]
  rest = data_files[1:]

  reference_data = pickle.load(open(reference, 'rb'))

  assert value_column in reference_data
  variance_column = None
  for column in ['%s.variance' % value_column,
                 value_column.replace('value', 'variance')]:
    if column in reference_data:
      variance_column = column
      break

  assert(variance_column)

  # construct XYZ pixel position search target

  reference_data = strip_not_integrated(reference_data)

  # keep only data with I/sig(I) > 3 for reference
  strong = (reference_data[value_column] > 3 * flex.sqrt(
    reference_data[variance_column]))

  print 'Keeping %d strong reflections of %d' % (strong.count(True),
                                                 len(reference_data))

  reference_data = reference_data.select(strong)

  xyz = reference_data['xyzcal.px'].as_double()
  ann = ann_adaptor(data=xyz, dim=3, k=1)

  for qpno, query_pickle in enumerate(rest):
    x = flex.double()
    y = flex.double()
    fout = open('matches%02d.dat' % qpno, 'w')
    query_data = strip_not_integrated(pickle.load(open(query_pickle, 'rb')))
    qxyz = query_data['xyzcal.px'].as_double()
    ann.query(qxyz)
    matches = 0
    for j, refl in enumerate(query_data):
      rrefl = reference_data[ann.nn[j]]
      if refl['miller_index'] == rrefl['miller_index']:
        fout.write('%d %d %d ' % refl['miller_index'] +
                   '%f %f ' % (rrefl[value_column], rrefl[variance_column]) +
                   '%f %f ' % (refl[value_column], refl[variance_column]) +
                   '%f %f %f\n' % refl['xyzcal.px'])
        matches += 1
        x.append(rrefl[value_column])
        y.append(refl[value_column])
    print 'For %s matched %d/%d' % (query_pickle, matches, len(query_data))
    fout.close()
    from matplotlib import pyplot
    pyplot.scatter(x.as_numpy_array(), y.as_numpy_array())
    pyplot.show()

if __name__ == '__main__':
  import sys
  value_column = sys.argv[1]
  saturation_analysis(sys.argv[2:], value_column)
