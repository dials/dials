# LIBTBX_SET_DISPATCHER_NAME dev.dials.semisynthetic_variance_analysis

from __future__ import division

def semisynthetic_variance_analysis(semisynthetic_integrated_data_files):
  import cPickle as pickle
  from logging import info
  from dials.array_family import flex
  from dials.util.add_hash import add_hash, dehash

  integrated_data_sets = [pickle.load(open(data_file, 'rb')) for
                          data_file in semisynthetic_integrated_data_files]

  # first prepare the data files i.e. remove partials, keep only integrated
  # reflections, add the hash column, add weight column

  hash_set = None

  for j, integrated_data in enumerate(integrated_data_sets):
    sel = integrated_data.get_flags(integrated_data.flags.integrated)
    integrated_data = integrated_data.select(sel)
    sel = integrated_data['partiality'] > 0.99
    integrated_data = integrated_data.select(sel)
    integrated_data = add_hash(integrated_data)
    integrated_data['intensity.prf.weight'] = 1.0 / integrated_data['intensity.prf.variance']
    integrated_data['intensity.sum.weight'] = 1.0 / integrated_data['intensity.sum.variance']
    print 'For file %s have %d reflections' % (semisynthetic_integrated_data_files[j], integrated_data.size())
    integrated_data = add_hash(integrated_data)
    hashes = set(integrated_data['hash'])
    if hash_set is None:
      hash_set = hashes
    else:
      hash_set = hash_set.intersection(hashes)

  print 'Union of all reflections has %d elements' % len(hash_set)

  # now figure a list of hashed reflections

if __name__ == '__main__':
  import sys
  semisynthetic_variance_analysis(sys.argv[1:])
