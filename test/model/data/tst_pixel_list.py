
class Test(object):

    def __init__(self):
        pass

    def run(self):
        self.tst_construct()
        self.tst_add_image()

    def tst_construct(self):

        from dials.model.data import PixelList
        from scitbx.array_family import flex
        from random import randint
        size = (2000, 2000)
        sf = 10
        pl = PixelList(size, sf)
        assert(pl.size() == size)
        assert(pl.first_frame() == sf)
        assert(pl.last_frame() == sf)
        assert(pl.num_frames() == 0)
        assert(pl.frame_range() == (sf, sf))

        frame_range = (10, 20)
        values = flex.int(range(100))
        coords = flex.vec3_int(100)
        for i in range(100):
            coords[i] = (
                randint(10, 20-1),
                randint(0, 2000-1),
                randint(0, 2000-1))

        pl = PixelList(size, frame_range, values, coords)
        assert(pl.size() == size)
        assert(pl.first_frame() == frame_range[0])
        assert(pl.last_frame() == frame_range[1])
        assert(pl.num_frames() == frame_range[1] - frame_range[0])
        assert(pl.frame_range() == frame_range)
        assert(len(pl.values()) == 100)
        assert(len(pl.coords()) == 100)
        for i in range(100):
            assert(pl.values()[i] == values[i])
            assert(pl.coords()[i] == coords[i])
        print 'OK'

    def tst_add_image(self):
        from dials.model.data import PixelList
        from scitbx.array_family import flex
        from random import randint
        size = (2000, 2000)
        sf = 10
        pl = PixelList(size, sf)

        count = 0
        for i in range(3):
            image = flex.random_int_gaussian_distribution(size[0]*size[1], 100, 5)
            mask = flex.random_bool(size[0]*size[1], 0.5)
            image.reshape(flex.grid(size))
            mask.reshape(flex.grid(size))
            count += len(mask.as_1d().select(mask.as_1d()))
            pl.add_image(image, mask)
        assert(len(pl.values()) == count)

if __name__ == '__main__':
    test = Test()
    test.run()
