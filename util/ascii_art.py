from __future__ import division

def spot_counts_per_image_plot(reflections, char='*', width=60, height=10, scan_range=None):
  import math
  from dials.array_family import flex

  assert isinstance(char, basestring)
  assert len(char) == 1

  x,y,z = reflections['xyzobs.px.value'].parts()

  # z-coordinates 0-1 lie on the first image. Internally this is image number 0,
  # but for the user this is image number 1. For the purposes of the histogram
  # all spots that are exactly on the first image (z=0.5) should be mapped to
  # image #1, therefore:
  z = z + 0.5

  if scan_range is None:
    max_z = int(math.ceil(flex.max(z)))
    min_z = int(math.floor(flex.min(z)))
  else:
    min_z, max_z = scan_range

  z_range = max_z - min_z
  if z_range == 1:
    return ''

  width = min(z_range, width)
  z_step = z_range / width

  counts = flex.double()

  sel = (z < min_z + z_step)
  counts.append(sel.count(True))
  for i in range(1, width-1):
    sel = ((z >= (min_z + i * z_step)) & (z < min_z + (i + 1) * z_step))
    counts.append(sel.count(True))
  sel = (z >= max_z - z_step)
  counts.append(sel.count(True))
# print list(z)
# for i in range(width):
#   print i, (min_z + i * z_step, min_z + (i + 1) * z_step), counts[i]

  max_count = flex.max(counts)
  total_counts = flex.sum(counts)
  counts *= (height/max_count)
  counts = counts.iround()

  rows = []
  rows.append('%i spots found on %i images (max %i / bin)' %(
    total_counts, z_range + 1, max_count))

  for i in range(height, 0, -1):
    row = []
    for j in range(len(counts)):
      if counts[j] > (i-1):
        row.append(char)
      else:
        row.append(' ')
    rows.append(''.join(row))

  padding = width - len(str(min_z)) - len(str(max_z))
  rows.append('%i%s%i' % (min_z,
    (' ' if padding < 7 else 'image').center(padding),
    max_z))
  return '\n'.join(rows)

if __name__ == '__main__':
  import sys
  from libtbx import easy_pickle

  for arg in sys.argv[1:]:
    print arg
    print spot_counts_per_image_plot(easy_pickle.load(arg))
    print
