

if __name__ == '__main__':
  from dials.array_family import flex
  import sys

  # Get the reflections from file
  table = flex.reflection_table.from_pickle(sys.argv[1])

  # Get the column data
  z = table['xyzcal.px'].parts()[2]
  I = table['intensity.sum.value']
  sig_I = flex.sqrt(table['intensity.sum.variance'])
  bg = table['background']

  # Calculate the average I, I / sig_I for each frame
  maxz = int(flex.max(z)) + 1
  count = [0] * maxz
  avr_I = [0] * maxz
  avr_IsigI = [0] * maxz
  avr_bg = [0] * maxz
  for zz, II, sig_II, bgg in zip(z, I, sig_I, bg):
    zz = int(zz)
    count[zz] += 1
    avr_I[zz] += II
    avr_IsigI[zz] += II / sig_II
    avr_bg[zz] += bgg

  for i in range(len(count)):
    avr_I[i] /= count[i]
    avr_IsigI[i] /= count[i]
    avr_bg[i] /= count[i]

  # Write out the averages to text file
  with open("average_intensity_per_frame.txt", "w") as outfile:
    for z, (II, sigII, bgg) in enumerate(zip(avr_I, avr_IsigI, avr_bg)):
      outfile.write('%d, %f, %f, %f\n' % (z, II, sigII, bgg))

  # Plot some graphs
  from matplotlib import pylab
  pylab.plot(avr_I)
  pylab.plot(avr_IsigI)
  pylab.plot(avr_bg)
  pylab.show()
