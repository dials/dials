import cPickle as pickle
from dials.array_family import flex
import math

#table = pickle.load(open('integrated_reflections_table.pickle', 'rb'))

table = pickle.load(open('../../../../../../xrd_2d_data/testing_detailed_tutorial_script_n3_20may_2015/indexed.pickle', 'rb'))
#show_reflections(table)

old_example = '''
print dir(table)
print len(table)
print "table.ncols() =", table.ncols()
print "table.nrows() =", table.nrows()
print "______________________________________________________________"
print "table.cols() =", table.cols()
print "table.rows() =", table.rows()
'''
print "table.keys() =", table.keys()
# Try iterating through table rows
for i in range(5):
  row = table[i]
  print row

  print "_______________________________"

example_from_test = '''
# Append some rows to the table
    row = { 'col1' : 10 }
    c1 = c1 + [10]
    c2 = c2 + [0]
    c3 = c3 + ['']
    table.append(row)
    assert(table.nrows() == 21)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'


'''


efisient_way = '''
its = table['intensity.cor.value']
for i in its:
  print i
#'''

stron_ref_table = flex.reflection_table()
for row in table.rows():
  h = row['miller_index']
  i_c = row['intensity.cor.value']
  i_r = row['intensity.sum.value']

  i_c_var = row['intensity.cor.variance']
  i_r_var = row['intensity.sum.variance']
  if( i_c > math.sqrt(i_c_var) and i_r > math.sqrt(i_r_var) ):
    stron_ref_table.append(row)
    #print h, i_c , i_r


code_to_reproduce = '''
  r_list = pickle.load(open(integrate_pkl, 'rb'))

  strong_reflections = []

  for r in r_list:
    if r.intensity > math.sqrt(r.intensity_variance):
      strong_reflections.append(r)
'''

not_in_use = '''
from dials.model.data import ReflectionList
rlist = ReflectionList.from_table(table)

for i in range(10):
  print rlist[i]
'''

