from dials.array_family import flex


'''
import numpy
data2d = numpy.zeros((3, 3), dtype = numpy.float64)

data2d[:, :] = 15
data2d[1:2, 1:2] = 50

for row in range(3):
  for col in range(3):
    data2d[row, col] += row * 2
    data2d[row, col] += col * 2

ref_table = flex.reflection_table()
'''


shoebox = flex.shoebox(5)
#panel = flex.size_t(1)
ref_table = flex.reflection_table()
shoebox.data = flex.double(flex.grid(6, 6, 6))
ref_table['shoebox'] = shoebox
#ref_table['panel'] = panel




#efisient_way = '''
its = ref_table['shoebox']
for arr in its:
  print arr.data.as_numpy_array()
#'''

example_from_other_code = '''
stron_ref_table = flex.reflection_table()
for row in table.rows():
  h = row['miller_index']
  i_c = row['intensity.cor.value']
  i_r = row['intensity.raw.value']

  i_c_var = row['intensity.cor.variance']
  i_r_var = row['intensity.raw.variance']
  if( i_c > math.sqrt(i_c_var) and i_r > math.sqrt(i_r_var) ):
    stron_ref_table.append(row)
'''

example_from_other_code = '''
    from random import randint
    from dials.model.data import Shoebox
    from scitbx.array_family import flex

    for i in range(1000):

      x0 = randint(0, 1000)
      y0 = randint(0, 1000)
      z0 = randint(0, 1000)
      x1 = randint(1, 10) + x0
      y1 = randint(1, 10) + y0
      z1 = randint(1, 10) + z0

      try:
        shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
        assert(shoebox.is_consistent() == False)
        shoebox.allocate()
        assert(shoebox.is_consistent() == True)
        shoebox.data = flex.double(flex.grid(20,20, 20))
        assert(shoebox.is_consistent() == False)
        shoebox.deallocate()
        assert(shoebox.is_consistent() == False)
      except Exception, e:
        print x0, y0, z0, x1, y1, z1
        raise

    # Test passed
    print 'OK'
'''


