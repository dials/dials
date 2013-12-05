
from dials.array_family import flex
from dials.framework.table import column_table
table = column_table()
table["col1"] = flex.double(10)
table["col2"] = flex.int(10)
col1 = table["col1"]
col2 = table["col2"]

print len(col1), len(col2)

col1.append(10)
col1.append(10)
col1.append(10)
col1.append(10)

print len(col1), len(col2)

col2.append(10)
col2.append(10)
col2.append(10)
col2.append(10)

print len(col1), len(col2)

col1.resize(10000000)

print len(col1), len(col2), len(table), table.nrows()

