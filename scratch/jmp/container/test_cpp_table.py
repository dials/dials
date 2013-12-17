

from dials.framework.table import column_table
from scitbx.array_family import flex

table = column_table()
table['c1'] = flex.int(10)
table['c2'] = flex.double(20)
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
