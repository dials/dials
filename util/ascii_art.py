from __future__ import division

def spot_counts_per_image_plot(reflections, char='*', width=60, height=10):
  import math
  from dials.array_family import flex

  x,y,z = reflections['xyzobs.px.value'].parts()
  max_z = int(math.ceil(flex.max(z)))
  min_z = int(math.floor(flex.min(z)))

  z_range = max_z - min_z
  width = min(z_range, width)

  z_step = z_range / width

  counts = flex.double()
  for i in range(width):
    sel = ((z >= (min_z + i * z_step)) & (z < min_z + (i + 1) * z_step))
    counts.append(sel.count(True))

  max_count = flex.max(counts)
  total_counts = flex.sum(counts)
  counts *= (height/max_count)
  counts = counts.iround()

  rows = []
  rows.append('%i spots found on %i images (max %i / bin)' %(
    total_counts, z_range, max_count))

  for i in range(10, 0, -1):
    row = []
    for j in range(len(counts)):
      if counts[j] > (i-1):
        row.append('*')
      else:
        row.append(' ')
    rows.append(''.join(row))

  from libtbx.math_utils import iceil
  first_image = '%i' %(min_z+1)
  last_image = '%i' %max_z
  if width < 10:
    padding = ' ' * width - len(first_image) - len(last_image)
  else:
    words = 'image'
    padding1 = ' ' * (
      iceil(width/2) - len(first_image) - iceil(len(words)/2))
    padding2 = ' ' * (
      iceil(width/2) - len(last_image) - iceil(len(words)/2))
    padding = padding1 + words + padding2

  while len(first_image + padding + last_image) < width:
    padding += ' '
  rows.append(first_image + padding + last_image)

  return '\n'.join(rows)


if __name__ == '__main__':
  import sys
  from libtbx import easy_pickle

  for arg in sys.argv[1:]:
    print arg
    print spot_counts_per_image_plot(easy_pickle.load(arg))
    print
