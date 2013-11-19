from __future__ import division

from table import Table

class ReflectionTable(Table):

  def spatial_index_map(self, column):
    pass


if __name__ == '__main__':
  from cctbx.array_family import flex

  # Create a reflection table with miller indices and intensities
  table = ReflectionTable([
    ('miller_index',
     flex.miller_index([(0, 0, 1), (0, 0, 2), (0, 0, 3), (0, 0, 2), (0, 0, 1)])
    ),
    ('intensity',
     flex.double([1, 2, 3, 4, 5])
    )
  ])

  print '\nGet the list of intensities for each miller index'
  indexer = table.index_map('miller_index')
  for key, indices in indexer.iteritems():
    print "Miller Index: {0}, Intensities: {1}".format(
      key, [table[i].intensity for i in indices])
