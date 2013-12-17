

from dials.framework.table import column_table
from scitbx.array_family import flex

table = column_table()
table['c1'] = flex.int(range(10))
table['c2'] = flex.double(range(20))
table['c3'] = flex.std_string(30)

print "Keys"
print table.keys()

print "Values"
print table.values()

print "Items"
print table.items()

print "Iterkeys"
for key in table.iterkeys():
  print key

print "IterValues"
for value in table.itervalues():
  print value

print "IterItems"
for key, value in table.iteritems():
  print key, value

print "IterRows"
for row in table.iterrows():
  print row

print "IndexRows"
print table[10]

print "Slice Table"
new_table = table[5:15]
print len(new_table)
for row in new_table:
  print row

print "Set Slice"
table[15:20] = new_table
for row in table:
  print row
