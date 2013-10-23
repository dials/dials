from __future__ import division

class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_allocate()
        self.tst_offset()
        self.tst_size()
        self.tst_consistent()
        self.tst_complete()


    def tst_allocate(self):
        from random import randint
        from dials.model.data import PartialShoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0 + 2
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11 + 2
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            shoebox.allocate()
            assert(shoebox.data.all() == (z11 - z00, y1 - y0, x1 - x0))
            shoebox.deallocate()
            assert(shoebox.data.all() == (0, 0, 0))

        # Test passed
        print 'OK'

    def tst_offset(self):
        from random import randint
        from dials.model.data import PartialShoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0 + 2
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11 + 2
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            assert(shoebox.xoffset() == x0)
            assert(shoebox.yoffset() == y0)
            assert(shoebox.zoffset() == z0)
            assert(shoebox.offset() == (z0, y0, x0))
            assert(shoebox.partial_zoffset() == z00)
            assert(shoebox.partial_offset() == (z00, y0, x0))

        # Test passed
        print 'OK'

    def tst_size(self):
        from random import randint
        from dials.model.data import PartialShoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0 + 2
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11 + 2
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            assert(shoebox.xsize() == x1 - x0)
            assert(shoebox.ysize() == y1 - y0)
            assert(shoebox.zsize() == z1 - z0)
            assert(shoebox.size() == (z1 - z0, y1 - y0, x1 - x0))
            assert(shoebox.partial_zsize() == z11 - z00)
            assert(shoebox.partial_size() == (z11 - z00, y1 - y0, x1 - x0))

        # Test passed
        print 'OK'

    def tst_consistent(self):
        from random import randint
        from dials.model.data import PartialShoebox
        from scitbx.array_family import flex

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0 + 2
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11 + 2
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            assert(shoebox.is_consistent() == False)
            shoebox.allocate()
            assert(shoebox.is_consistent() == True)
            shoebox.data = flex.int(flex.grid(10, 10, 10))
            assert(shoebox.is_consistent() == False)
            shoebox.deallocate()
            assert(shoebox.is_consistent() == False)

        # Test passed
        print 'OK'

    def tst_complete(self):
        from random import randint
        from dials.model.data import PartialShoebox
        from scitbx.array_family import flex

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0 + 2
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11 + 2
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            assert(shoebox.is_complete() == False)

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            z00 = z0
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z11 = randint(1, 10) + z00
            z1 = z11
            bbox = (x0, x1, y0, y1, z0, z1)
            zrange = z00, z11

            shoebox = PartialShoebox(bbox, zrange)
            assert(shoebox.is_complete() == True)

        # Test passed
        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()
