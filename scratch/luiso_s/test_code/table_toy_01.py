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
  n += 1
  img = arr.data
  img = flex.double(flex.grid(3, 3, 3))

  for row in range(3):
    for col in range(3):
      for fra in range(3):
        img[row, col, fra] = row + col + fra + n * 9
  arr.data = img

its = ref_table['intensity.raw.value']
i_var = ref_table['intensity.raw.variance']

for i in range(num_ref):
  its[i] = (i + 1) * 11
  i_var[i] = (i + 1) * 12

print ">>>>>>>>>>>>>>>>>>>>>>>>>    printing before integrrating         <<<<<<<<"

iterate = ref_table['shoebox']
for arr in iterate:
  np_img = arr.data.as_numpy_array()
  print np_img
  print ">>"

iterate = ref_table['intensity.raw.value']
for n_its in iterate:
  print n_its
print ">>>"
iterate = ref_table['intensity.raw.variance']
for n_i_v in iterate:
  print n_i_v
