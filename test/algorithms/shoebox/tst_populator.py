from __future__ import division

class Test(object):

    def __init__(self):
        ''' Create the maps and reflection list'''
        from scitbx.array_family import flex

        # Create the mask, gainmap and dark map
        xsize, ysize = (2000, 2000)
        self.mask = flex.bool(flex.grid(ysize, xsize), True)
        self.gain_map = flex.double(flex.grid(ysize, xsize), 1.0)
        self.dark_map = flex.double(flex.grid(ysize, xsize), 0.0)

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
        self.tst_initialize()
        self.tst_add_image()
        self.tst_image_mask()

    def tst_initialize(self):
        ''' Test that the indices of each pixel are the same'''
        from dials.algorithms import shoebox
        from collections import defaultdict
        from scitbx.array_family import flex

        # Create the populator
        populator = shoebox.Populator(self.reflections, self.mask,
                                      self.gain_map, self.dark_map)

        # Create a dictionary of flex arrays
        frames_to_reflection = defaultdict(flex.int)

        # For each reflection, Find the frames which it spans and copy an
        # index into the frame -> reflection list
        for i, r in enumerate(self.reflections):
            if r.is_valid():
                f0 = r.bounding_box[4]
                f1 = r.bounding_box[5]
                for f in range(f0, f1):
                    frames_to_reflection[f].append(i)

        # Ensure both give the same indices
        for k, v in frames_to_reflection.iteritems():
            index1 = v
            index2 = populator.indices(k)
            assert(index1 == index2)


        # Test passed
        print 'OK'

    def tst_add_image(self):
        ''' Add a load of images and then check that the values of the
        shoebox pixels are correct.'''
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Create the populator
        shoebox.allocate(self.reflections)
        populator = shoebox.Populator(self.reflections, self.mask,
                                      self.gain_map, self.dark_map)

        # Add all the images to the shoeboxes
        images = []
        for i in range(10):
            image = flex.random_int_gaussian_distribution(2000*2000, 100, 10)
            image.reshape(flex.grid(2000, 2000))
            images.append(image)
            populator.add_image(image, i)

        # The accuracy
        eps = 1e-7

        # Check all reflection shoeboxes have the correct values
        for r in self.reflections:
            x0, x1, y0, y1, z0, z1 = r.bounding_box
            x00, x11, y00, y11, z00, z11 = r.bounding_box
            profile = r.shoebox

            # Clip coordinates
            if x0 < 0: x0 = 0
            if y0 < 0: y0 = 0
            if z0 < 0: z0 = 0
            if x1 > 2000: x1 = 2000
            if y1 > 2000: y1 = 2000
            if z1 > 10: z1 = 10

            # Ensure that the values are either all zero or all the
            # same as the image at the same coordinates
            if x1 <= x0 or y1 <= y0 or z1 <= z0:
                assert(abs(profile).all_le(eps))
            else:
                for k in range(z0, z1):
                    for j in range(y0, y1):
                        for i in range(x0, x1):
                            kk = k - z00
                            jj = j - y00
                            ii = i - x00
                            value1 = images[k][j, i]
                            value2 = profile[kk, jj, ii]
                            assert(abs(value1 - value2) <= eps)

        # Test passed
        print 'OK'

    def tst_image_mask(self):
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Create the populator and deallocate shoeboxes
        populator = shoebox.Populator(self.reflections, self.mask,
                                      self.gain_map, self.dark_map)

        # Create the mask and ensure that no pixels are set that shouldn't be
        # and that all pixels are set that should be
        for i in range(10):
            mask1 = populator.image_mask(i, (0, 0))
            assert(mask1.all()[0] == 2000 and mask1.all()[1] == 2000)

            mask2 = flex.bool(flex.grid(2000, 2000), False)
            for r in self.reflections:
                x0, x1, y0, y1, z0, z1 = r.bounding_box
                if z0 <= i and z1 > i:
                    if y0 < 0: y0 = 0
                    if x0 < 0: x0 = 0
                    if y1 > 2000: y1 = 2000
                    if x1 > 2000: x1 = 2000
                    if x1 > x0 and y1 > y0:
                        mask2[y0:y1,x0:x1] = flex.bool(
                            flex.grid(y1-y0,x1-x0), True)
            assert(mask1 == mask2)

        # Test passed
        print 'OK'

if __name__ == '__main__':
    test = Test()
    test.run()
