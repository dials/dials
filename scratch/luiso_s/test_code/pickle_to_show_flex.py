import cPickle as pickle
from dials.array_family import flex
import math

table = pickle.load(open('../../../../../../xrd_2d_data/testing_detailed_tutorial_script_n3_20may_2015/indexed.pickle', 'rb'))
#show_reflections(table)


print "table.keys() =", table.keys()
# Try iterating through table rows
for i in range(5):
  row = table[i]
  print row

  print "_______________________________"


from dials.viewer.slice_viewer import show_3d
'''
flex_dat_frst_lst = []
flex_dat_seg_lst = []

for nm in [98296, 102718, 100823, 9313, 106348, 31151, 81228, 71186, 59599]:

  flex_dat_frst_lst.append(table[nm]['shoebox'].data)
  flex_dat_seg_lst.append(table[nm]['shoebox'].mask)

'''
flex_dat_frst_lst = []
flex_dat_seg_lst = []

for nm in range(59540, 59750):

  flex_dat_frst_lst.append(table[nm]['shoebox'].data)
  flex_dat_seg_lst.append(table[nm]['shoebox'].mask)

'''
#to see data from a single shoebox
show_3d(flex_dat_frst_lst[row])

#to see data and mask from a single shoebox
show_3d(flex_dat_frst_lst[row], flex_dat_seg_lst[0])
'''
#to see data from a set shoeboxes
show_3d(flex_dat_frst_lst)

#to see data and mask from a set shoeboxes
show_3d(flex_dat_frst_lst, flex_dat_seg_lst)


for row in table.rows():
  shoebox_data = row['shoebox'].data
  shoebox_mask = row['shoebox'].mask

