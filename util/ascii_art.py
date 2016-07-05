from __future__ import division

def spot_counts_per_image_plot(reflections, char='*', width=60, height=10):
  from dials.array_family import flex

  if len(reflections) == 0:
    return '\n'

  assert isinstance(char, basestring)
  assert len(char) == 1

  x,y,z = reflections['xyzobs.px.value'].parts()

  min_z = flex.min(z)
  max_z = flex.max(z)

  # image numbers to display on x-axis label
  xlab = (int(round(min_z + 0.5)), int(round(max_z + 0.5)))
  # estimate the total number of images
  image_count = xlab[1] - xlab[0] + 1

  z_range = max_z - min_z + 1
  if z_range <= 1:
    return '%i spots found on 1 image' %len(reflections)

  width = int(min(z_range, width))
  z_step = z_range / width
  z_bound = min_z + z_step - 0.5
# print [round(i * 10) / 10 for i in sorted(z)]

  counts = flex.double()

  sel = (z < z_bound)
  counts.append(sel.count(True))
# print 0, ('-', z_bound), sel.count(True)
  for i in range(1, width-1):
    sel = ((z >= z_bound) & (z < (z_bound + z_step)))
    counts.append(sel.count(True))
#   print i, (z_bound, z_bound + z_step), sel.count(True)
    z_bound += z_step
  sel = (z >= z_bound)
# print i + 1, (z_bound, '-'), sel.count(True)
  counts.append(sel.count(True))

  max_count = flex.max(counts)
  total_counts = flex.sum(counts)
  assert total_counts == len(z)
  counts *= (height/max_count)
  counts = counts.iround()

  rows = []
  rows.append('%i spots found on %i images (max %i / bin)' %(
    total_counts, image_count, max_count))

  for i in range(height, 0, -1):
    row = []
    for j, c in enumerate(counts):
      if c > (i - 1):
        row.append(char)
      else:
        row.append(' ')
    rows.append(''.join(row))

  padding = width - len(str(xlab[0])) - len(str(xlab[1]))
  rows.append('%i%s%i' % (xlab[0],
    (' ' if padding < 7 else 'image').center(padding),
    xlab[1]))
  return '\n'.join(rows)

if __name__ == '__main__':
  import sys
  from libtbx import easy_pickle

  for arg in sys.argv[1:]:
    print arg
    print spot_counts_per_image_plot(easy_pickle.load(arg))
    print
