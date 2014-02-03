
def assert_almost(a, b, eps=1e-7):
  if isinstance(a, tuple):
    assert(len(a) == len(b))
    for aa, bb in zip(a, b):
      assert_almost(aa, bb)
  else:
    assert(abs(a - b) < eps)

class Test:
  def __init__(self):
    pass

  def run(self):
    self.tst_to_table()
    self.tst_from_table()

  def generate_reflection_list(self, n):
    from dials.model.data import ReflectionList
    import random
    rand_int = lambda: random.randint(0, 100)
    rand_double = lambda: random.uniform(0, 100)
    rand_bool = lambda: random.choice([True, False])
    rand_vec3_int = lambda: (rand_int(), rand_int(), rand_int())
    rand_vec3_double = lambda: (rand_double(), rand_double(), rand_double())
    rand_int6 = lambda: rand_vec3_int() + rand_vec3_int()
    rand_vec2_double = lambda: (rand_double(), rand_double())

    rlist = ReflectionList(n)
    for r in rlist:
      r.miller_index = rand_vec3_int()
      r.status = rand_int()
      r.crystal = rand_int()
      r.panel = rand_int()
      r.entering = rand_bool()
      r.bounding_box = rand_int6()
      r.beam_vector = rand_vec3_double()
      r.centroid_position = rand_vec3_double()
      r.centroid_variance = rand_vec3_double()
      r.intensity = rand_double()
      r.intensity_variance = rand_double()
      r.corrected_intensity = rand_double()
      r.corrected_intensity_variance = rand_double()
      r.image_coord_mm = rand_vec2_double()
      r.image_coord_px = rand_vec2_double()
      r.rotation_angle = rand_double()
      r.frame_number = rand_double()

    return rlist

  def generate_reflection_table(self, n):
    from dials.array_family import flex
    import random
    rand_int = lambda: random.randint(0, 100)
    rand_double = lambda: random.uniform(0, 100)
    rand_bool = lambda: random.choice([True, False])
    rand_vec3_int = lambda: (rand_int(), rand_int(), rand_int())
    rand_vec3_double = lambda: (rand_double(), rand_double(), rand_double())
    rand_int6 = lambda: rand_vec3_int() + rand_vec3_int()
    rand_vec2_double = lambda: (rand_double(), rand_double())

    flex_rand_int = lambda n: flex.int([rand_int() for i in range(n)])
    flex_rand_double = lambda n:flex.double([rand_double() for i in range(n)])
    flex_rand_bool = lambda n:flex.bool([rand_bool() for i in range(n)])
    flex_rand_vec3_int = lambda n:flex.miller_index([rand_vec3_int() for i in range(n)])
    flex_rand_vec3_double = lambda n:flex.vec3_double([rand_vec3_double() for i in range(n)])
    flex_rand_int6 = lambda n:flex.int6([rand_int6() for i in range(n)])
    flex_rand_vec2_double = lambda n:flex.vec2_double([rand_vec2_double() for i in range(n)])

    return flex.reflection_table([
      ('hkl', flex_rand_vec3_int(n)),
      ('flags', flex_rand_int(n)),
      ('id', flex_rand_int(n)),
      ('panel', flex_rand_int(n)),
      ('entering', flex_rand_bool(n)),
      ('shoebox.bbox', flex_rand_int6(n)),
      ('s1', flex_rand_vec3_double(n)),
      ('xyzcal.px', flex_rand_vec3_double(n)),
      ('xyzcal.mm', flex_rand_vec3_double(n)),
      ('xyzobs.px.value', flex_rand_vec3_double(n)),
      ('xyzobs.px.variance', flex_rand_vec3_double(n)),
      ('intensity.raw.value', flex_rand_double(n)),
      ('intensity.raw.variance', flex_rand_double(n)),
      ('intensity.cor.value', flex_rand_double(n)),
      ('intensity.cor.variance', flex_rand_double(n))
    ])

  def check_values(self, t, l):

    assert(t.nrows() == len(l))
    for r1, r2 in zip(t.rows(), l):
      assert(r1['hkl'] == r2.miller_index)
      assert(r1['flags'] == r2.status)
      assert(r1['id'] == r2.crystal)
      assert(r1['panel'] == r2.panel_number)
      assert(r1['entering'] == r2.entering)
      assert(r1['shoebox.bbox'] == r2.bounding_box)
      assert_almost(r1['s1'], r2.beam_vector)
      assert_almost(r1['xyzobs.px.value'], r2.centroid_position)
      assert_almost(r1['xyzobs.px.variance'], r2.centroid_variance)
      assert_almost(r1['intensity.raw.value'], r2.intensity)
      assert_almost(r1['intensity.raw.variance'], r2.intensity_variance)
      assert_almost(r1['intensity.cor.value'], r2.corrected_intensity)
      assert_almost(r1['intensity.cor.variance'], r2.corrected_intensity_variance)

  def tst_to_table(self):

    # Create a reflection list
    n = 100
    rlist = self.generate_reflection_list(n)
    table = rlist.to_table()
    assert(table.is_consistent())
    assert(table.nrows() == n)
    assert(table.ncols() == 15)
    self.check_values(table, rlist)
    print 'OK'

  def tst_from_table(self):
    from dials.model.data import ReflectionList
    n = 100
    table = self.generate_reflection_table(n)
    rlist = ReflectionList.from_table(table)
    self.check_values(table, rlist)
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
