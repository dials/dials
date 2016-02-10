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

    self.tst_init()
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
    self.tst_flags()
    self.tst_copy()
    self.tst_extract_shoeboxes()
    self.tst_split_by_experiment_id()
    self.tst_split_indices_by_experiment_id()
    self.tst_split_partials()
    self.tst_split_partials_with_shoebox()
    self.tst_find_overlapping()

  def tst_init(self):
    from dials.array_family import flex

    # test default
    table = flex.reflection_table()
    assert(table.is_consistent())
    assert(table.nrows() == 0)
    assert(table.ncols() == 0)
    assert(table.empty())
    print 'Ok'

    # test with nrows
    table = flex.reflection_table(10)
    assert(table.is_consistent())
    assert(table.nrows() == 10)
    assert(table.ncols() == 0)
    assert(table.empty())
    print 'OK'

    # test with valid columns
    table = flex.reflection_table([
      ('col1', flex.int(10)),
      ('col2', flex.double(10)),
      ('col3', flex.std_string(10))])
    assert(table.is_consistent())
    assert(table.nrows() == 10)
    assert(table.ncols() == 3)
    assert(not table.empty())
    print 'OK'

    # test with invalid columns
    try:
      table = flex.reflection_table([
        ('col1', flex.int(10)),
        ('col2', flex.double(20)),
        ('col3', flex.std_string(10))])
      assert(false)
    except Exception:
      pass
    print 'OK'

  def tst_resizing(self):
    from dials.array_family import flex

    # Create a table with 2 empty columns
    table = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table0 = flex.reflection_table()
    table1 = flex.reflection_table()
    table2 = flex.reflection_table()
    table1['col1'] = flex.int(c1)
    table1['col2'] = flex.double(c2)
    table2['col3'] = flex.std_string(c3)

    # Update from zero columns
    table0.update(table1)
    assert(table0.is_consistent())
    assert(table0.nrows() == 10)
    assert(table0.ncols() == 2)
    print 'OK'

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
    table2 = flex.reflection_table()
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
    table = flex.reflection_table()
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
    table1 = flex.reflection_table()
    table2 = flex.reflection_table()
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
    table1 = flex.reflection_table()
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
    table2 = flex.reflection_table()
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
    table2 = flex.reflection_table()
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

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table1 = flex.reflection_table()
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
    table1 = flex.reflection_table()
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
    table1 = flex.reflection_table()
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
    table1 = flex.reflection_table()
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

  def tst_flags(self):

    from dials.array_family import flex

    # Create a table with flags all 0
    table = flex.reflection_table()
    table['flags'] = flex.size_t(5, 0)

    # Get all the flags
    f1 = table.get_flags(table.flags.predicted)
    assert(f1.count(True) == 0)
    print 'OK'

    # Set some flags
    mask = flex.bool([True, True, False, False, True])
    table.set_flags(mask, table.flags.predicted)
    f1 = table.get_flags(table.flags.predicted)
    assert(f1.count(True) == 3)
    assert(all(f11 == f22 for f11, f22 in zip(f1, mask)))
    f2 = table.get_flags(table.flags.predicted | table.flags.observed)
    assert(f2.count(True) == 0)
    print 'OK'

    # Unset the flags
    mask = flex.bool(5, True)
    table.unset_flags(mask, table.flags.predicted | table.flags.observed)
    f1 = table.get_flags(table.flags.predicted)
    assert(f1.count(True) == 0)
    flags = table['flags']
    assert(all(f == 0 for f in flags))
    print 'OK'

    # Set multiple flags
    mask = flex.bool([True, True, False, False, True])
    table.set_flags(mask, table.flags.predicted | table.flags.observed)
    f1 = table.get_flags(table.flags.predicted)
    f2 = table.get_flags(table.flags.observed)
    assert(f1.count(True) == 3)
    assert(f2.count(True) == 3)
    mask = flex.bool([False, True, True, True, False])
    table.set_flags(mask, table.flags.integrated)
    f1 = table.get_flags(table.flags.predicted)
    f2 = table.get_flags(table.flags.observed)
    f3 = table.get_flags(table.flags.integrated)
    f4 = table.get_flags(table.flags.integrated | table.flags.predicted)
    assert(f1.count(True) == 3)
    assert(f2.count(True) == 3)
    assert(f3.count(True) == 3)
    assert(f4.count(True) == 1)
    print 'OK'

    # Get where any are set
    f1 = table.get_flags(table.flags.predicted, all=False)
    f2 = table.get_flags(table.flags.observed, all=False)
    f3 = table.get_flags(table.flags.integrated, all=False)
    f4 = table.get_flags(table.flags.integrated | table.flags.predicted, all=False)
    assert(f1.count(True) == 3)
    assert(f2.count(True) == 3)
    assert(f3.count(True) == 3)
    assert(f4.count(True) == 5)
    print 'OK'

  def tst_serialize(self):

    from dials.array_family import flex

    # The columns as lists
    c1 = list(range(10))
    c2 = list(range(10))
    c3 = ['a', 'b', 'c', 'd', 'e', 'f', 'g', 'i', 'j', 'k']

    # Create a table with some elements
    table = flex.reflection_table()
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

  def tst_copy(self):
    import copy
    from dials.array_family import flex

    # Create a table
    table = flex.reflection_table([
      ('col1', flex.int(range(10)))])

    # Make a shallow copy of the table
    shallow = copy.copy(table)
    shallow['col2'] = flex.double(range(10))
    assert(table.ncols() == 2)
    assert(table.is_consistent())
    print 'OK'

    # Make a deep copy of the table
    deep = copy.deepcopy(table)
    deep['col3'] = flex.std_string(10)
    assert(table.ncols() == 2)
    assert(deep.ncols() == 3)
    assert(table.is_consistent())
    assert(deep.is_consistent())

    table2 = table.copy()
    table2['col3'] = flex.std_string(10)
    assert(table.ncols() == 2)
    assert(table2.ncols() == 3)
    assert(table.is_consistent())
    assert(table2.is_consistent())
    print 'OK'

  def tst_extract_shoeboxes(self):
    from dials.array_family import flex
    from random import randint, seed
    from dials.algorithms.shoebox import MaskCode
    import sys
    seed(0)

    reflections = flex.reflection_table()
    reflections['panel'] = flex.size_t()
    reflections['bbox'] = flex.int6()

    npanels = 2
    width = 1000
    height = 1000
    nrefl = 10000
    frame0 = 10
    frame1 = 100
    nrefl = 1000

    for i in range(nrefl):
      xs = randint(5, 10)
      ys = randint(5, 10)
      x0 = randint(-xs+1, width-1)
      y0 = randint(-ys+1, height-1)
      z0 = randint(frame0, frame1-1)
      x1 = x0 + xs
      y1 = y0 + ys
      z1 = min([z0 + randint(1, 10), frame1])
      assert(x1 > x0)
      assert(y1 > y0)
      assert(z1 > z0)
      assert(z0 >= frame0 and z1 <= frame1)
      bbox = (x0, x1, y0, y1, z0, z1)
      reflections.append({
        "panel" : randint(0,1),
        "bbox" : bbox,
      })

    reflections['shoebox'] = flex.shoebox(
      reflections['panel'],
      reflections['bbox'])
    reflections['shoebox'].allocate()

    class FakeImageSet(object):
      def __init__(self):
        from dials.array_family import flex
        self.data = flex.int(range(height*width))
        self.data.reshape(flex.grid(height,width))
      def get_array_range(self):
        return (frame0, frame1)
      def get_detector(self):
        class FakeDetector(object):
          def __len__(self):
            return npanels
          def __getitem__(self, index):
            class FakePanel(object):
              def get_trusted_range(self):
                return (-1, 1000000)
            return FakePanel()
        return FakeDetector()
      def __len__(self):
        return frame1 - frame0
      def __getitem__(self, index):
        f = frame0+index
        return (self.data + f*1, self.data + f*2)
      def get_corrected_data(self, index):
        f = frame0+index
        return (self.data + f*1, self.data + f*2)
      def get_mask(self, index):
        image = self.get_corrected_data(index)
        return tuple(im >= 0 for im in image)
    imageset = FakeImageSet()

    stdout = sys.stdout
    class DevNull(object):
      def write(self, *args):
        pass
      def flush(self):
        pass
    sys.stdout = DevNull()
    reflections.extract_shoeboxes(imageset)
    sys.stdout = stdout

    for i in range(len(reflections)):
      sbox = reflections[i]["shoebox"]
      assert(sbox.is_consistent())
      mask = sbox.mask
      data = sbox.data
      bbox = sbox.bbox
      panel = sbox.panel
      x0, x1, y0, y1, z0, z1 = bbox
      for z in range(z1 - z0):
        for y in range(y1 - y0):
          for x in range(x1 - x0):
            v1 = data[z,y,x]
            m1 = mask[z,y,x]
            if (x0 + x >= 0 and y0 + y >= 0 and
                x0 + x < width and y0 + y < height):
              v2 = imageset.data[y+y0,x+x0] + (z+z0)*(panel+1)
              m2 = MaskCode.Valid
              assert(v1 == v2)
              assert(m1 == m2)
            else:
              assert(v1 == 0)
              assert(m1 == 0)

    print 'OK'

  def tst_split_by_experiment_id(self):
    from dials.array_family import flex
    r = flex.reflection_table()
    r['id'] = flex.int()
    for i in range(100):
      r.append({"id" : 0})
      r.append({"id" : 1})
      r.append({"id" : 2})
      r.append({"id" : 3})
      r.append({"id" : 5})
    result = r.split_by_experiment_id()
    assert(len(result) == 5)
    for res, exp in zip(result, [0, 1, 2, 3, 5]):
      assert(len(res) == 100)
      assert(res['id'].count(exp) == 100)
    print 'OK'

  def tst_split_indices_by_experiment_id(self):
    from dials.array_family import flex
    r = flex.reflection_table()
    r['id'] = flex.int()
    for i in range(100):
      r.append({"id" : 0})
      r.append({"id" : 1})
      r.append({"id" : 2})
      r.append({"id" : 3})
      r.append({"id" : 5})
    index_list = r.split_indices_by_experiment_id(6)
    assert(len(index_list) == 6)
    for index, exp, num in zip(
        index_list,
        [0, 1, 2, 3, 4, 5],
        [100, 100, 100, 100, 0, 100]):
      assert(len(index) == num)
      assert(r.select(index)['id'].count(exp) == num)
    print 'OK'

  def tst_split_partials(self):
    from dials.array_family import flex
    from random import randint, uniform
    r = flex.reflection_table()
    r['value1'] = flex.double()
    r['value2'] = flex.int()
    r['value3'] = flex.double()
    r['bbox'] = flex.int6()
    expected = []
    for i in range(100):
      x0 = randint(0, 100)
      x1 = x0 + randint(1, 10)
      y0 = randint(0, 100)
      y1 = y0 + randint(1, 10)
      z0 = randint(0, 100)
      z1 = z0 + randint(1, 10)
      v1 = uniform(0, 100)
      v2 = randint(0, 100)
      v3 = uniform(0, 100)
      r.append({
        'value1' : v1,
        'value2' : v2,
        'value3' : v3,
        'bbox' : (x0, x1, y0, y1, z0, z1)
      })
      for z in range(z0, z1):
        expected.append({
          'value1' : v1,
          'value2' : v2,
          'value3' : v3,
          'bbox' : (x0, x1, y0, y1, z, z+1),
          'partial_id' : i,
        })

    r.split_partials()
    assert(len(r) == len(expected))
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
      assert(abs(r1['value1'] - r2['value1']) < EPS)
      assert(r1['value2'] == r2['value2'])
      assert(abs(r1['value3'] - r2['value3']) < EPS)
      assert(r1['bbox'] == r2['bbox'])
      assert(r1['partial_id'] == r2['partial_id'])

    print 'OK'

  def tst_split_partials_with_shoebox(self):
    from dials.array_family import flex
    from random import randint, uniform
    from dials.model.data import Shoebox
    r = flex.reflection_table()
    r['value1'] = flex.double()
    r['value2'] = flex.int()
    r['value3'] = flex.double()
    r['bbox'] = flex.int6()
    r['panel'] = flex.size_t()
    r['shoebox'] = flex.shoebox()
    expected = []
    for i in range(100):
      x0 = randint(0, 100)
      x1 = x0 + randint(1, 10)
      y0 = randint(0, 100)
      y1 = y0 + randint(1, 10)
      z0 = randint(0, 100)
      z1 = z0 + randint(1, 10)
      v1 = uniform(0, 100)
      v2 = randint(0, 100)
      v3 = uniform(0, 100)
      sbox = Shoebox(0, (x0, x1, y0, y1, z0, z1))
      sbox.allocate()
      assert(sbox.is_consistent())
      w = x1 - x0
      h = y1 - y0
      for z in range(z0, z1):
        for y in range(y0, y1):
          for x in range(x0, x1):
            sbox.data[z-z0,y-y0,x-x0] = x+y*w+z*w*h
      r.append({
        'value1' : v1,
        'value2' : v2,
        'value3' : v3,
        'bbox' : (x0, x1, y0, y1, z0, z1),
        'panel' : 0,
        'shoebox' : sbox
      })
      for z in range(z0, z1):
        sbox = Shoebox(0, (x0, x1, y0, y1, z, z+1))
        sbox.allocate()
        assert(sbox.is_consistent())
        w = x1 - x0
        h = y1 - y0
        for y in range(y0, y1):
          for x in range(x0, x1):
            sbox.data[0,y-y0,x-x0] = x+y*w+z*w*h
        expected.append({
          'value1' : v1,
          'value2' : v2,
          'value3' : v3,
          'bbox' : (x0, x1, y0, y1, z, z+1),
          'partial_id' : i,
          'panel' : 0,
          'shoebox' : sbox
        })

    r.split_partials_with_shoebox()
    assert(len(r) == len(expected))
    EPS = 1e-7
    for r1, r2 in zip(r, expected):
      assert(abs(r1['value1'] - r2['value1']) < EPS)
      assert(r1['value2'] == r2['value2'])
      assert(abs(r1['value3'] - r2['value3']) < EPS)
      assert(r1['bbox'] == r2['bbox'])
      assert(r1['partial_id'] == r2['partial_id'])
      assert(r1['panel'] == r2['panel'])
      assert(r1['shoebox'].data.as_double().as_1d().all_approx_equal(
        r2['shoebox'].data.as_double().as_1d()))

    print 'OK'

  def tst_find_overlapping(self):
    from dials.array_family import flex
    from random import randint, uniform
    from dials.model.data import Shoebox
    N = 10000
    r = flex.reflection_table(N)
    r['bbox'] = flex.int6(N)
    r['panel'] = flex.size_t(N)
    r['id'] = flex.int(N)
    r['imageset_id'] = flex.int(N)
    for i in range(N):
      x0 = randint(0, 100)
      x1 = randint(1, 10) + x0
      y0 = randint(0, 100)
      y1 = randint(1, 10) + y0
      z0 = randint(0, 100)
      z1 = randint(1, 10) + z0
      panel = randint(0,2)
      pid = randint(0,2)
      r['bbox'][i] = (x0,x1,y0,y1,z0,z1)
      r['panel'][i] = panel
      r['id'][i] = pid
      r['imageset_id'][i] = pid

    def is_overlap(b0, b1, border):
      b0 = b0[0]-border,b0[1]+border,b0[2]-border,b0[3]+border,b0[4]-border,b0[5]+border
      b1 = b1[0]-border,b1[1]+border,b1[2]-border,b1[3]+border,b1[4]-border,b1[5]+border
      if not (b1[0] > b0[1] or
              b1[1] < b0[0] or
              b1[2] > b0[3] or
              b1[3] < b0[2] or
              b1[4] > b0[5] or
              b1[5] < b0[4]):
        return True
      return False

    for i in [0, 2, 5]:
      overlaps = r.find_overlaps(border=i)
      for item in overlaps.edges():
        i0 = overlaps.source(item)
        i1 = overlaps.target(item)
        r0 = r[i0]
        r1 = r[i1]
        p0 = r0['panel']
        p1 = r1['panel']
        b0 = r0['bbox']
        b1 = r1['bbox']
        j0 = r0['imageset_id']
        j1 = r1['imageset_id']
        assert j0 == j1
        assert p0 == p1
        assert is_overlap(b0,b1,i)



    print 'OK'


if __name__ == '__main__':
  from dials.test import cd_auto
  with cd_auto(__file__):
    test = Test()
    test.run()
