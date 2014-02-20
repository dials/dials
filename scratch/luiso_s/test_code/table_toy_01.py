from dials.array_family import flex
num_ref = 3
ref_table = flex.reflection_table()

shoebox = flex.shoebox(num_ref)
ref_table['shoebox'] = shoebox

intensity = flex.double(num_ref)
ref_table['intensity.raw.value'] = intensity

intensity_var = flex.double(num_ref)
ref_table['intensity.raw.variance'] = intensity_var

iterate = ref_table['shoebox']
n = 0
for arr in iterate:
  img = flex.double(flex.grid(3, 3, 3))
  bkg = flex.double(flex.grid(3, 3, 3))
  msk = flex.int(flex.grid(3, 3, 3))
  for row in range(3):
    for col in range(3):
      for fra in range(3):
        img[row, col, fra] = row + col + fra + n * 9
        bkg[row, col, fra] = 0.0
        msk[row, col, fra] = 3
  n += 1
  msk[1, 1, 1] = 5
  tmp_i = n * n * n * 3
  img[1, 1, 1] += tmp_i
  print "intensity must be =", tmp_i
  arr.data = img[:, :, :]
  arr.background = bkg[:, :, :]
  arr.mask = msk[:, :, :]

its = ref_table['intensity.raw.value']
i_var = ref_table['intensity.raw.variance']

for i in range(num_ref):
  its[i] = (i + 1) * 11
  i_var[i] = (i + 1) * 12

print ">>>>>>>>>>>>>>>>>>>>>>>>>    printing before integrating         <<<<<<<<"

iterate = ref_table['shoebox']
for arr in iterate:
  np_img = arr.data.as_numpy_array()
  print np_img
<<<<<<< Updated upstream
  np_img = arr.background.as_numpy_array()
  print np_img
  np_img = arr.mask.as_numpy_array()
  print np_img

  print ">>"

iterate = ref_table['intensity.raw.value']
for n_its in iterate:
  print n_its
print ">>>"
iterate = ref_table['intensity.raw.variance']
for n_i_v in iterate:
  print n_i_v


'''
from dials.algorithms.background.curved_background_subtractor \
 import layering_and_background_modl
layering_and_background_modl(ref_table)


from dials.algorithms.background.inclined_background_subtractor \
 import layering_and_background_plane
layering_and_background_plane(ref_table)
'''

from dials.algorithms.background.flat_background_subtractor \
 import layering_and_background_avg
layering_and_background_avg(ref_table)


print ">>>>>>>>>>>>>>>>>>>>>>>>>    printing during integrating         <<<<<<<<"
print ">>>>>>>>>>>>>>>>>>>>>>>>>      printing with background           <<<<<<<<"
iterate = ref_table['shoebox']
for arr in iterate:
  np_img = arr.data.as_numpy_array()
  print np_img
  np_img = arr.background.as_numpy_array()
  print np_img
  np_img = arr.mask.as_numpy_array()
  print np_img
  print ">>"

iterate = ref_table['intensity.raw.value']
for n_its in iterate:
  print n_its
print ">>>"
iterate = ref_table['intensity.raw.variance']
for n_i_v in iterate:
  print n_i_v
=======
>>>>>>> Stashed changes
