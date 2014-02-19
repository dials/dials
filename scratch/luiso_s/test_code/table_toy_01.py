from dials.array_family import flex


shoebox = flex.shoebox(5)

intensity = 6.7

ref_table = flex.reflection_table()
ref_table['shoebox'] = shoebox
#ref_table['intensity'] = intensity

iterate = ref_table['shoebox']
for arr in iterate:
  img = arr.data
  img = flex.double(flex.grid(3, 3, 3))

  for row in range(3):
    for col in range(3):
      for fra in range(3):
        img[row, col, fra] = row + col + fra
  print img.as_numpy_array()
  arr.data = img

iterate = ref_table['shoebox']
for arr in iterate:

  np_img = arr.data.as_numpy_array()
  print np_img




