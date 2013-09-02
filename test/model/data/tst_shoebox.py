
class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_allocate()
        self.tst_offset()
        self.tst_size()
        self.tst_consistent()
        self.tst_is_bbox_within_image_volume()
        self.tst_does_bbox_contain_bad_pixels()
        self.tst_count_mask_values()

    def tst_allocate(self):
        from random import randint
        from dials.model.data import Shoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            shoebox.allocate()
            assert(shoebox.data.all() == (z1 - z0, y1 - y0, x1 - x0))
            assert(shoebox.mask.all() == (z1 - z0, y1 - y0, x1 - x0))
            shoebox.deallocate()
            assert(shoebox.data.all() == (0, 0, 0))
            assert(shoebox.mask.all() == (0, 0, 0))

        # Test passed
        print 'OK'

    def tst_offset(self):
        from random import randint
        from dials.model.data import Shoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            assert(shoebox.xoffset() == x0)
            assert(shoebox.yoffset() == y0)
            assert(shoebox.zoffset() == z0)
            assert(shoebox.offset() == (z0, y0, x0))

        # Test passed
        print 'OK'

    def tst_size(self):
        from random import randint
        from dials.model.data import Shoebox

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            assert(shoebox.xsize() == x1 - x0)
            assert(shoebox.ysize() == y1 - y0)
            assert(shoebox.zsize() == z1 - z0)
            assert(shoebox.size() == (z1 - z0, y1 - y0, x1 - x0))

        # Test passed
        print 'OK'

    def tst_consistent(self):
        from random import randint
        from dials.model.data import Shoebox
        from scitbx.array_family import flex

        for i in range(10):

            x0 = randint(0, 1000)
            y0 = randint(0, 1000)
            z0 = randint(0, 1000)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            assert(shoebox.is_consistent() == False)
            shoebox.allocate()
            assert(shoebox.is_consistent() == True)
            shoebox.data.resize(flex.grid(10, 10, 10))
            assert(shoebox.is_consistent() == False)
            shoebox.deallocate()
            assert(shoebox.is_consistent() == False)

        # Test passed
        print 'OK'

    def tst_is_bbox_within_image_volume(self):

        from dials.model.data import Shoebox

        isize = (1000, 1000)
        srange = (0, 100)

        shoebox = Shoebox((10, 20, 10, 20, 10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == True)
        shoebox = Shoebox((-10, 20, 10, 20, 10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)
        shoebox = Shoebox((10, 20, -10, 20, 10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)
        shoebox = Shoebox((10, 20, 10, 20, -10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)
        shoebox = Shoebox((10, 1020, 10, 20, 10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)
        shoebox = Shoebox((10, 20, 10, 1020, 10, 20))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)
        shoebox = Shoebox((10, 20, 10, 20, 10, 1020))
        assert(shoebox.is_bbox_within_image_volume(isize, srange) == False)

        # Test passed
        print 'OK'

    def tst_does_bbox_contain_bad_pixels(self):
        from scitbx.array_family import flex
        from dials.model.data import Shoebox
        from random import randint

        mask = flex.bool(flex.grid(100, 100), True)
        for j in range(100):
            for i in range(40, 60):
                mask[j,i] = False
                mask[i,j] = False

        for i in range(1000):
            x0 = randint(0, 90)
            y0 = randint(0, 90)
            z0 = randint(0, 90)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))

            res1 = shoebox.does_bbox_contain_bad_pixels(mask)
            res2 = False
            if x0 >= 40 and x0 < 60:
                res2 = True
            if x1 > 40 and x1 <= 60:
                res2 = True
            if y0 >= 40 and y0 < 60:
                res2 = True
            if y1 > 40 and y1 <= 60:
                res2 = True

            assert(res1 == res2)

        # Test passed
        print 'OK'

    def tst_count_mask_values(self):

        from dials.model.data import Shoebox
        from random import randint, sample

        for i in range(10):
            x0 = randint(0, 90)
            y0 = randint(0, 90)
            z0 = randint(0, 90)
            x1 = randint(1, 10) + x0
            y1 = randint(1, 10) + y0
            z1 = randint(1, 10) + z0

            shoebox = Shoebox((x0, x1, y0, y1, z0, z1))
            shoebox.allocate()
            maxnum = len(shoebox.mask)
            num = randint(1, maxnum)
            indices = sample(list(range(maxnum)), num)
            value = (1 << 2)
            for i in indices:
                shoebox.mask[i] = value

            assert(shoebox.count_mask_values(value) == num)

        # Test passed
        print 'OK'


if __name__ == '__main__':
    test = Test()
    test.run()
