from __future__ import division

def assert_exception(func):
  try:
    func()
    assert(False)
  except Exception:
    pass

class Test(object):
  def __init__(self):
    pass

  def run(self):

    self.tst_resizing()
    self.tst_iteration()
    self.tst_row_operations()
    self.tst_slicing()
    self.tst_updating()
    self.tst_select()
    self.tst_set_selected()
    self.tst_serialize()
    self.tst_delete()
    self.tst_del_selected()

  def tst_resizing(self):
    from dials.array_family import flex

    # Create a table with 2 empty columns
    table = ReflectionTable()
    assert(table.empty())
    table['col1'] = flex.int()
    table['col2'] = flex.double()
    assert(table.nrows() == 0)
    assert(table.ncols() == 2)
    assert(not table.empty())
    assert('col1' in table)
    assert('col2' in table)
    assert('col3' not in table)
    print 'OK'

    # Create a table with 2 columns and 10 rows
    table = ReflectionTable()
    table['col1'] = flex.int(10)
    table['col2'] = flex.double(10)
    assert(table.nrows() == 10)
    assert(table.ncols() == 2)
    print 'OK'

    # Add an extra column with the wrong size (throw)
    try:
      table['col3'] = flex.std_string(20)
      assert(False)
    except Exception:
      pass
    assert(table.nrows() == 10)
    assert(table.ncols() == 2)
    assert(table.is_consistent())
    assert(len(table['col1']) == 10)
    assert(len(table['col2']) == 10)
    print 'OK'

    # Resize the table (should resize all columns)
    table.resize(50)
    assert(table.nrows() == 50)
    assert(table.ncols() == 2)
    assert(table.is_consistent())
    assert(len(table['col1']) == 50)
    assert(len(table['col2']) == 50)
    print 'OK'

    # Make the table inconsistent
    table['col1'].resize(40)
    assert(not table.is_consistent())
    assert_exception(lambda: table.nrows())
    assert_exception(lambda: table.ncols())
    print 'OK'

    # Clear the table
    table.clear()
    assert(table.is_consistent())
    assert(table.empty())
    assert(table.nrows() == 0)
    assert(table.ncols() == 0)
    print 'OK'

  def tst_delete(self):
    from dials.array_family import flex

    # Test del item
    table = ReflectionTable()
    table['col1'] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    table['col2'] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    table['col3'] = flex.int([0, 1, 2, 3, 4, 5, 6, 7, 8, 9])
    del table['col3']
    assert(table.is_consistent())
    assert(table.nrows() == 10)
    assert(table.ncols() == 2)
    assert(not "col3" in table)
    print 'OK'

    # Test del row
    del table[5]
    assert(table.is_consistent())
    assert(table.nrows() == 9)
    assert(table.ncols() == 2)
    assert(all(a==b for a, b in zip(list(table['col1']),
      [0, 1, 2, 3, 4, 6, 7, 8, 9])))
    print 'OK'

    # Test del slice
    del table[0:10:2]
    assert(table.is_consistent())
    assert(table.nrows() == 4)
    assert(table.ncols() == 2)
    assert(all(a==b for a, b in zip(list(table['col1']),
      [1, 3, 6, 8])))
    print 'OK'

    # Test del slice
    del table[:]
    assert(table.is_consistent())
    assert(table.nrows() == 0)
    assert(table.ncols() == 2)
    print 'OK'

  def tst_row_operations(self):
    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = ReflectionTable()
    table['col1'] = flex.int(c1)
    table['col2'] = flex.double(c2)
    table['col3'] = flex.std_string(c3)

    # Extend the table
    table.extend(table)
    c1 = c1 * 2
    c2 = c2 * 2
    c3 = c3 * 2
    assert(table.nrows() == 20)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

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

    row = { 'col2' : 11 }
    c1 = c1 + [0]
    c2 = c2 + [11]
    c3 = c3 + ['']
    table.append(row)
    assert(table.nrows() == 22)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

    row = { 'col1' : 12, 'col2' : 12, 'col3' : 'l' }
    c1 = c1 + [12]
    c2 = c2 + [12]
    c3 = c3 + ['l']
    table.append(row)
    assert(table.nrows() == 23)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

    # Try inserting some rows
    row = { 'col1' : -1 }
    c1.insert(5, -1)
    c2.insert(5, 0)
    c3.insert(5, '')
    table.insert(5, row)
    assert(table.nrows() == 24)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

    row = { 'col1' : -2, 'col2' : -3, 'col3' : 'abc' }
    c1.insert(2, -2)
    c2.insert(2, -3)
    c3.insert(2, 'abc')
    table.insert(2, row)
    assert(table.nrows() == 25)
    assert(table.ncols() == 3)
    assert(table.is_consistent())
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

    # Try iterating through table rows
    for i in range(table.nrows()):
      row = table[i]
      assert(row['col1'] == c1[i])
      assert(row['col2'] == c2[i])
      assert(row['col3'] == c3[i])
    print 'OK'

    # Trying setting some rows
    row = { 'col1' : 100 }
    table[2] = row
    assert(table[2]['col1'] == 100)
    assert(table[2]['col2'] == c2[2])
    assert(table[2]['col3'] == c3[2])

    row = { 'col1' : 1000, 'col2' : 2000, 'col3' : 'hello' }
    table[10] = row
    assert(table[10]['col1'] == 1000)
    assert(table[10]['col2'] == 2000)
    assert(table[10]['col3'] == 'hello')
    print 'OK'

  def tst_iteration(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = ReflectionTable()
    table['col1'] = flex.int(c1)
    table['col2'] = flex.double(c2)
    table['col3'] = flex.std_string(c3)

    # Try iterating keys
    k = []
    for key in table.keys():
      k.append(key)
    assert(len(k) == 3)
    assert(k.count('col1') == 1)
    assert(k.count('col2') == 1)
    assert(k.count('col3') == 1)
    print 'OK'

    # Try iterating columns
    k = []
    c = []
    for key, col in table.cols():
      k.append(key)
      c.append(col)
    assert(len(k) == 3)
    assert(k.count('col1') == 1)
    assert(k.count('col2') == 1)
    assert(k.count('col3') == 1)
    print 'OK'

    # Try iterating rows
    for row1, row2 in zip(table.rows(), zip(c1, c2, c3)):
      assert(row1['col1'] == row2[0])
      assert(row1['col2'] == row2[1])
      assert(row1['col3'] == row2[2])
    print 'OK'

  def tst_slicing(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = ReflectionTable()
    table['col1'] = flex.int(c1)
    table['col2'] = flex.double(c2)
    table['col3'] = flex.std_string(c3)

    # Try forward slicing
    new_table = table[2:7:2]
    assert(new_table.ncols() == 3)
    assert(new_table.nrows() == 3)
    assert(new_table.is_consistent())
    c11 = c1[2:7:2]
    c22 = c2[2:7:2]
    c33 = c3[2:7:2]
    assert(all(a == b for a, b in zip(new_table['col1'], c11)))
    assert(all(a == b for a, b in zip(new_table['col2'], c22)))
    assert(all(a == b for a, b in zip(new_table['col3'], c33)))
    print 'OK'

    # Try backward slicing
    new_table = table[7:2:-2]
    assert(new_table.ncols() == 3)
    assert(new_table.nrows() == 3)
    assert(new_table.is_consistent())
    c11 = c1[7:2:-2]
    c22 = c2[7:2:-2]
    c33 = c3[7:2:-2]
    assert(all(a == b for a, b in zip(new_table['col1'], c11)))
    assert(all(a == b for a, b in zip(new_table['col2'], c22)))
    assert(all(a == b for a, b in zip(new_table['col3'], c33)))
    print 'OK'

    # Try setting forward slicing
    table[2:7:2] = new_table
    assert(table.ncols() == 3)
    assert(table.nrows() == 10)
    assert(table.is_consistent())
    c1[2:7:2] = c11
    c2[2:7:2] = c22
    c3[2:7:2] = c33
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

    # Try setting backward slicing
    table[7:2:-2] = new_table
    assert(table.ncols() == 3)
    assert(table.nrows() == 10)
    assert(table.is_consistent())
    c1[7:2:-2] = c11
    c2[7:2:-2] = c22
    c3[7:2:-2] = c33
    assert(all(a == b for a, b in zip(table['col1'], c1)))
    assert(all(a == b for a, b in zip(table['col2'], c2)))
    assert(all(a == b for a, b in zip(table['col3'], c3)))
    print 'OK'

  def tst_updating(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table1 = ReflectionTable()
    table2 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table2['col3'] = flex.std_string(c3)

    # Update table1 with table2 columns
    table1.update(table2)
    assert(table1.is_consistent())
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 3)
    assert(table2.is_consistent())
    assert(table2.nrows() == 10)
    assert(table2.ncols() == 1)
    print 'OK'

    # Update trable1 with invalid table
    c3 = ['a', 'b', 'c']

    # Create a table with some elements
    table2 = ReflectionTable()
    table2['col3'] = flex.std_string(c3)
    try:
      table1.update(table2)
      assert(False)
    except Exception:
      pass

    assert(table1.is_consistent())
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 3)
    assert(table2.is_consistent())
    assert(table2.nrows() == 3)
    assert(table2.ncols() == 1)
    print 'OK'

  def tst_select(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = ReflectionTable()
    table['col1'] = flex.int(c1)
    table['col2'] = flex.double(c2)
    table['col3'] = flex.std_string(c3)

    # Select some columns
    new_table = table.select(('col1', 'col2'))
    assert(new_table.nrows() == 10)
    assert(new_table.ncols() == 2)
    assert(all(a == b for a, b in zip(new_table['col1'], c1)))
    assert(all(a == b for a, b in zip(new_table['col2'], c2)))
    print 'OK'

    # Select some columns
    new_table = table.select(flex.std_string(['col1', 'col2']))
    assert(new_table.nrows() == 10)
    assert(new_table.ncols() == 2)
    assert(all(a == b for a, b in zip(new_table['col1'], c1)))
    assert(all(a == b for a, b in zip(new_table['col2'], c2)))
    print 'OK'

    # Select some rows
    index = flex.size_t([0, 1, 5, 8, 9])
    cc1 = [c1[i] for i in index]
    cc2 = [c2[i] for i in index]
    cc3 = [c3[i] for i in index]
    new_table = table.select(index)
    assert(new_table.nrows() == 5)
    assert(new_table.ncols() == 3)
    assert(all(a == b for a, b in zip(new_table['col1'], cc1)))
    assert(all(a == b for a, b in zip(new_table['col2'], cc2)))
    assert(all(a == b for a, b in zip(new_table['col3'], cc3)))
    print 'OK'

    # Select some rows
    index = flex.bool([True, True, False, False, False,
                       True, False, False, True, True])
    new_table = table.select(index)
    assert(new_table.nrows() == 5)
    assert(new_table.ncols() == 3)
    assert(all(a == b for a, b in zip(new_table['col1'], cc1)))
    assert(all(a == b for a, b in zip(new_table['col2'], cc2)))
    assert(all(a == b for a, b in zip(new_table['col3'], cc3)))
    print 'OK'

  def tst_set_selected(self):

    from dials.array_family import flex
    from copy import deepcopy

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table1 = ReflectionTable()
    table2 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table2['col2'] = flex.double(c2)
    table2['col3'] = flex.std_string(c3)

    # Set selected columns
    table1.set_selected(('col3', 'col2'), table2)
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 3)
    assert(all(a == b for a, b in zip(table1['col1'], c1)))
    assert(all(a == b for a, b in zip(table1['col2'], c2)))
    assert(all(a == b for a, b in zip(table1['col3'], c3)))
    print 'OK'

    # Set selected columns
    table1 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1.set_selected(flex.std_string(['col3', 'col2']), table2)
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 3)
    assert(all(a == b for a, b in zip(table1['col1'], c1)))
    assert(all(a == b for a, b in zip(table1['col2'], c2)))
    assert(all(a == b for a, b in zip(table1['col3'], c3)))
    print 'OK'

    cc1 = list(range(10, 15))
    cc2 = list(range(10, 15))
    cc3 = ['l', 'm', 'n', 'o', 'p']

    # Set selected rows
    table2 = ReflectionTable()
    table2['col1'] = flex.int(cc1)
    table2['col2'] = flex.double(cc2)
    table2['col3'] = flex.std_string(cc3)

    index = flex.size_t([0, 1, 5, 8, 9])
    ccc1 = deepcopy(c1)
    ccc2 = deepcopy(c2)
    ccc3 = deepcopy(c3)
    for j, i in enumerate(index):
      ccc1[i] = cc1[j]
      ccc2[i] = cc2[j]
      ccc3[i] = cc3[j]
    table1.set_selected(index, table2)
    assert(all(a == b for a, b in zip(table1['col1'], ccc1)))
    assert(all(a == b for a, b in zip(table1['col2'], ccc2)))
    assert(all(a == b for a, b in zip(table1['col3'], ccc3)))
    print 'OK'

    # Set selected rows
    table2 = ReflectionTable()
    table2['col1'] = flex.int(cc1)
    table2['col2'] = flex.double(cc2)
    table2['col3'] = flex.std_string(cc3)

    flags = flex.bool([True, True, False, False, False,
                       True, False, False, True, True])
    table1.set_selected(index, table2)
    assert(all(a == b for a, b in zip(table1['col1'], ccc1)))
    assert(all(a == b for a, b in zip(table1['col2'], ccc2)))
    assert(all(a == b for a, b in zip(table1['col3'], ccc3)))
    print 'OK'

  def tst_del_selected(self):

    from dials.array_family import flex
    from copy import deepcopy

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table1 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table1['col3'] = flex.std_string(c3)

    # Del selected columns
    table1.del_selected(('col3', 'col2'))
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 1)
    assert("col1" in table1)
    assert("col2" not in table1)
    assert("col3" not in table1)
    assert(all(a == b for a, b in zip(table1['col1'], c1)))
    print 'OK'

    # Del selected columns
    table1 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table1['col3'] = flex.std_string(c3)
    table1.del_selected(flex.std_string(['col3', 'col2']))
    assert(table1.nrows() == 10)
    assert(table1.ncols() == 1)
    assert("col1" in table1)
    assert("col2" not in table1)
    assert("col3" not in table1)
    assert(all(a == b for a, b in zip(table1['col1'], c1)))
    print 'OK'

    # Del selected rows
    table1 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table1['col3'] = flex.std_string(c3)

    index = flex.size_t([0, 1, 5, 8, 9])
    index2 = range(10)
    for i in index:
      index2.remove(i)
    ccc1 = [c1[i] for i in index2]
    ccc2 = [c2[i] for i in index2]
    ccc3 = [c3[i] for i in index2]
    table1.del_selected(index)
    assert(table1.nrows() == len(ccc1))
    assert(all(a == b for a, b in zip(table1['col1'], ccc1)))
    assert(all(a == b for a, b in zip(table1['col2'], ccc2)))
    assert(all(a == b for a, b in zip(table1['col3'], ccc3)))
    print 'OK'

    # Del selected rows
    table1 = ReflectionTable()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table1['col3'] = flex.std_string(c3)

    flags = flex.bool([True, True, False, False, False,
                       True, False, False, True, True])
    table1.del_selected(index)
    assert(table1.nrows() == len(ccc1))
    assert(all(a == b for a, b in zip(table1['col1'], ccc1)))
    assert(all(a == b for a, b in zip(table1['col2'], ccc2)))
    assert(all(a == b for a, b in zip(table1['col3'], ccc3)))
    print 'OK'

  def tst_serialize(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = ReflectionTable()
    table['col1'] = flex.int(c1)
    table['col2'] = flex.double(c2)
    table['col3'] = flex.std_string(c3)

    # Pickle, then unpickle
    import cPickle as pickle
    obj = pickle.dumps(table)
    new_table = pickle.loads(obj)
    assert(new_table.is_consistent())
    assert(new_table.nrows() == 10)
    assert(new_table.ncols() == 3)
    assert(all(a == b for a, b in zip(new_table['col1'], c1)))
    assert(all(a == b for a, b in zip(new_table['col2'], c2)))
    assert(all(a == b for a, b in zip(new_table['col3'], c3)))
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
