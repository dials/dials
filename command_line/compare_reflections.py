from __future__ import absolute_import, division, print_function

# LIBTBX_SET_DISPATCHER_NAME dev.dials.compare_reflections

class CompareReflections(object):
  ''' A class to compare reflections. '''

  def __init__(self, refl1, refl2):
    ''' Initialise the class. '''
    self.refl1 = refl1
    self.refl2 = refl2

  def intensities(self):
    ''' Compare the intensities. '''
    from dials.array_family import flex

    # Sort by resolution
    d = self.refl1['d']
    index = flex.size_t(reversed(sorted(range(len(d)), key=lambda x: d[x])))
    self.refl1.reorder(index)
    self.refl2.reorder(index)

    # Get the intensities
    I1 = self.refl1['intensity.sum.value']
    I2 = self.refl2['intensity.sum.value']
    S1 = flex.sqrt(self.refl1['intensity.sum.variance'])
    S2 = flex.sqrt(self.refl2['intensity.sum.variance'])
    xyz1 = self.refl1['xyzcal.px']
    xyz2 = self.refl2['xyzcal.px']

    # Compute chunked statistics
    corr = []
    R = []
    scale = []
    res = []
    for i in range(len(self.refl1) // 1000):

      # Get the chunks of data
      a = i * 1000
      b = (i+1) * 1000
      II1 = I1[a:b]
      II2 = I2[a:b]
      res.append(d[a])

      # Compute the mean and standard deviation per chunk
      mv1 = flex.mean_and_variance(II1)
      mv2 = flex.mean_and_variance(II2)
      m1 = mv1.mean()
      m2 = mv2.mean()
      s1 = mv1.unweighted_sample_standard_deviation()
      s2 = mv2.unweighted_sample_standard_deviation()

      # compute the correlation coefficient
      r = (1/(len(II1) - 1))*sum(((II1[j] - m1) / s1) * ((II2[j] - m2) / s2)
          for j in range(len(II1)))
      corr.append(r)

      # Compute the scale between the chunks
      s = sum(II1) / sum(II2)
      scale.append(s)

      # Compute R between the chunks
      r = sum(abs(abs(II1[j]) - abs(s * II2[j])) for j in range(len(II1))) \
        / sum(abs(II1[j]) for j in range(len(II1)))
      R.append(r)

    from matplotlib import pylab
    pylab.plot(corr, label="CC")
    pylab.plot(R, label="R")
    pylab.plot(scale, label="K")
    pylab.legend()
    pylab.show()


if __name__ == '__main__':

  from optparse import OptionParser
  from dials.util.command_line import Command
  from dials.array_family import flex
  from annlib_ext import AnnAdaptor as ann_adaptor
  import libtbx.load_env

  usage = "usage: %s [options] reflections1.pickle reflections2.pickle" \
    % libtbx.env.dispatcher_name
  parser = OptionParser(usage)

  # Parse the command line arguments
  (options, args) = parser.parse_args()

  # Ensure there are only two reflection lists
  if len(args) != 2:
    parser.print_help()
    exit(0)

  # Read the first batch of reflections
  Command.start('Reading reflections from %s' % args[0])
  refl1 = flex.reflection_table.from_pickle(args[0])
  mask = flex.bool(xyz == (0, 0, 0) for xyz in refl1['xyzobs.px.value'])
  refl1.del_selected(mask)
  Command.end('Read %d reflections from %s' % (len(refl1), args[0]))

  # Read the second batch of reflections
  Command.start('Reading reflections from %s' % args[1])
  refl2 = flex.reflection_table.from_pickle(args[1])
  mask = refl2['intensity.sum.value'] <= 0.0
  refl2.del_selected(mask)
  mask = refl2['intensity.sum.value']**2 < refl2['intensity.sum.variance']
  refl2.del_selected(mask)
  Command.end('Read %d reflections from %s' % (len(refl2), args[1]))

  # perform the match
  Command.start('Find matching reflections')
  hkl1 = refl1['miller_index']
  hkl2 = refl2['miller_index']
  xyz1 = refl1['xyzcal.px']
  xyz2 = refl2['xyzcal.px']

  # Do the nn match
  ann = ann_adaptor(xyz1.as_double().as_1d(), 3)
  ann.query(xyz2.as_double().as_1d())

  # Select only those with matching hkl
  index = flex.size_t(i for i in ann.nn)
  hkl11 = hkl1.select(index)
  flags = hkl11 == hkl2
  index = index.select(flags)
  refl1 = refl1.select(index)
  refl2 = refl2.select(flags)
  Command.end('Found %d matching reflections' % len(refl1))

  # Do the comparison
  compare = CompareReflections(refl1, refl2)
  compare.intensities()
