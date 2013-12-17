

from dials.model.data import reflection_table
from scitbx.array_family import flex
from time import time
table = reflection_table()
table['c1'] = flex.int()
table['c2'] = flex.double()
table.resize(1000000)

st = time()

keys, cols = zip(*list(table.cols()))
for c in zip(*cols):
  d = dict(zip(keys, c))
print time() - st

st = time()
for row in table.rows():
  pass
print time() - st

#table['c3'] = flex.std_string(30)

#print "Keys"
#print table.keys()

#print "Items"
#print table.cols()

#print "IterRows"
#for row in table.rows():
#  print row

#print "IndexRows"
#print table[10]

#print "Slice Table"
#new_table = table[5:15]
#print len(new_table)
#for row in new_table:
#  print row

#print "Set Slice"
#table[15:20] = new_table
#for row in table:
#  print row

#print "c1" in table
#print "c" in table

#table[10] = { 'c1' : 100 }
#table[11] = { 'c1' : 200, 'c3' : "Hello World" }

#print table[10]
#print table[11]

#table = column_table([
#  ("column_1", flex.int()),
#  ("column_2", flex.std_string())])

#print list(table.keys())

#table.append({ 'column_1' : 200, 'column_2' : "Hello World 1" })
#table.append({ 'column_1' : 300, 'column_2' : "Hello World 2" })
#table.append({ 'column_1' : 400, 'column_2' : "Hello World 3" })
#table.append({ 'column_1' : 500, 'column_2' : "Hello World 4" })

#for row in table.rows():
#  print row

#table.insert(2, { 'column_1' : 1000 })

#for row in table.rows():
#  print row

#print "Extend"
#table.extend(table)


#for row in table.rows():
#  print row

#new_table = column_table([
#  ("column_2", flex.int()),
#  ("column_3", flex.int()),
#  ("column_4", flex.std_string())])

#table.update(new_table)

#for row in table.rows():
#  print row


#print "Reorder"
#table = column_table()
#table["c1"] = flex.int([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#table["c2"] = flex.int([1, 2, 3, 4, 5, 6, 7, 8, 9, 10])
#table["c3"] = flex.int([9, 8, 7, 6, 5, 4, 3, 2, 1, 0])

#index = flex.size_t([9, 8, 7, 6, 4, 4, 3, 2, 1, 0])
#table.reorder(index)

#for row in table.rows():
#  print row

#print "Sort"
#table.sort("c1")

#for row in table.rows():
#  print row

#print "Sort"
#table.sort("c1", reverse=True)

#for row in table.rows():
#  print row

#print "Types"
#print table.types()
