from __future__ import division

class Test(object):

  def __init__(self):
    ''' Create the maps and reflection list'''
    from scitbx.array_family import flex

    # Generate some reflections
    self.reflections = self.generate_reflections(1000)

  def generate_reflections(self, num):
    ''' Generate some random reflections.'''
    from dials.model.data import ReflectionList
    from random import randint
    rlist = ReflectionList(num)
    for i in range(num):
      x0 = randint(-10, 2000)
      x1 = x0 + randint(1, 10)
      y0 = randint(-10, 2000)
      y1 = y0 + randint(1, 10)
      z0 = randint(-5, 10)
      z1 = z0 + randint(1, 10)
      rlist[i].bounding_box = x0, x1, y0, y1, z0, z1
    return rlist

  def run(self):
    ''' Run the tests'''
    self.tst_deallocate()
    self.tst_allocate()

  def tst_deallocate(self):
    ''' test the delallocate function deletes all the shoeboxes'''
    from dials.algorithms import shoebox

    # Allocate and deallocate the shoeboxes
    shoebox.allocate(self.reflections)
    shoebox.deallocate(self.reflections)

    # Check all reflection shoeboxes are empty
    for r in self.reflections:
      assert(len(r.shoebox) == 0)
      assert(len(r.shoebox_mask) == 0)
      assert(len(r.shoebox_background) == 0)

    # Test passed
    print 'OK'

  def tst_allocate(self):
    ''' Test that shoeboxes are all the right size.'''
    from dials.algorithms import shoebox

    # Allocate the shoeboxes
    shoebox.allocate(self.reflections)

    # Check all reflection shoeboxes are the right size
    for r in self.reflections:
      x0, x1, y0, y1, z0, z1 = r.bounding_box
      zs1 = z1 - z0
      ys1 = y1 - y0
      xs1 = x1 - x0
      zs2, ys2, xs2 = r.shoebox.all()
      zs3, ys3, xs3 = r.shoebox_mask.all()
      zs4, ys4, xs4 = r.shoebox_background.all()
      assert(zs1 == zs2 and zs1 == zs3 and zs1 == zs4)
      assert(ys1 == ys2 and ys1 == ys3 and ys1 == ys4)
      assert(xs1 == xs2 and xs1 == xs3 and xs1 == xs4)

    # Test passed
    print 'OK'

if __name__ == '__main__':
  test = Test()
  test.run()
