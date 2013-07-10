

class TestReciprocalSpaceTransformE3Fraction(object):

    def __init__(self, filename):

        from math import pi
        from dials.model.serialize import load
        from dials.algorithms.integration import \
            ReciprocalSpaceTransformE3Fraction
        from dials.algorithms.integration import BBoxCalculator

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        # Set the delta_divergence/mosaicity
        self.n_sigma = 3
        self.sigma_divergence = 0.060 * pi / 180
        self.mosaicity = 0.154 * pi / 180
        self.delta_divergence = self.n_sigma * self.sigma_divergence
        self.delta_mosaicity = self.n_sigma * self.mosaicity

        # Set the grid size
        self.grid_size = (4, 4, 4)

        # Create the E3 fraction object
        self.transform = ReciprocalSpaceTransformE3Fraction(
            self.scan,
            self.mosaicity,
            self.n_sigma,
            self.grid_size[2])

        # Create the bounding box calculator
        self.calculate_bbox = BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            self.delta_divergence,
            self.delta_mosaicity)

    def __call__(self):

        from scitbx import matrix
        from random import uniform
        from dials.algorithms.integration import FromBeamVectorToXds
        from dials.algorithms.integration import XdsCoordinateSystem
        from scitbx.array_family import flex

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        for i in range(1000):

            # Get random x, y, z
            x = uniform(0, 2000)
            y = uniform(0, 2000)
            z = uniform(-1000, 1000)

            # Get random s1, phi, panel
            s1 = matrix.col(self.detector.get_pixel_lab_coord(
                (x, y))).normalize() * s0_length
            phi = self.scan.get_angle_from_array_index(z, deg=False)
            panel = 0

            # Calculate the bounding box
            bbox = self.calculate_bbox(s1, phi, panel)
            x1, x2 = bbox[0], bbox[1]
            y1, y2 = bbox[2], bbox[3]
            z1, z2 = bbox[4], bbox[5]

            # Create the XDS coordinate system
            xcs = XdsCoordinateSystem(s0, s1, m2, phi)

            # Calculate the transform fraction
            fraction = self.transform(bbox[4:], phi, xcs.zeta)

            # Ensure the minimum and maximum are 0 < 1
            fmax = flex.max(fraction)
            fmin = flex.min(fraction)
            assert(fmax <= 1.0)
            assert(fmin >= 0.0)

            # Ensure the fraction for each image frame adds up to 1.0 for
            # all those frames completely within the grid
            for j in range(1, fraction.all()[0]-1):
                tot = flex.sum(fraction[j:j+1,:])
                assert(abs(tot - 1.0) < 1e-7)

            # Ensure the frames follow a progression through the grid. I.e,
            # check that values increase then decrease and don't jump around
            for j in range(fraction.all()[0]):
                f = fraction[j:j+1,:]
                last = f[0]
                rev = False
                for i in range(1,len(f)):
                  curr = f[1]
                  if rev == False:
                      if curr < last:
                          rev = True
                  else:
                      assert(curr <= last)
                  last = curr

        # Test passed
        print 'OK'


class TestReciprocalSpaceTransformDetectorLabCoords(object):


    def __init__(self, filename):
        from math import pi
        from dials.model.serialize import load
        from dials.algorithms.integration import \
            ReciprocalSpaceTransformDetectorLabCoords

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        # Set the number of sub-division
        self.n_div = 2

        # Create the calculator
        self.detector_s1 = ReciprocalSpaceTransformDetectorLabCoords()

    def __call__(self):

        from scitbx import matrix
        from random import randint

        # The detector beam vectors
        ds1 = self.detector_s1(self.detector, self.beam, self.n_div)

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        for k in range(1000):
            j = randint(0, ds1.all()[0])
            i = randint(0, ds1.all()[1])
            y = float(j + 0.5) / self.n_div
            x = float(i + 0.5) / self.n_div
            xyz = self.detector.get_pixel_lab_coord((x, y))
            s11 = matrix.col(xyz).normalize() * s0_length
            s12 = matrix.col(ds1[j,i])
            assert((s11 - s12).length() < 1e-7)

        print 'OK'


class TestReciprocalSpaceTransform(object):

    def __init__(self, filename):
        from math import pi
        from dials.model.serialize import load
        from dials.algorithms.integration import ReciprocalSpaceTransform
        from dials.algorithms.integration import BBoxCalculator

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        # Set the delta_divergence/mosaicity
        self.n_sigma = 5
        self.sigma_divergence = 0.060 * pi / 180
        self.mosaicity = 0.154 * pi / 180
        self.delta_divergence = self.n_sigma * self.sigma_divergence
        self.delta_mosaicity = self.n_sigma * self.mosaicity

        # Set the grid size
        self.grid_half_size = 5

        # Set the number of sub-division
        self.n_div = 5

        # Create the bounding box calculator
        self.calculate_bbox = BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            self.delta_divergence,
            self.delta_mosaicity)

        # Create the reciprocal space transform
        self.transform = ReciprocalSpaceTransform(
            self.beam, self.detector, self.gonio, self.scan,
            self.mosaicity, float(self.n_sigma),
            self.grid_half_size, self.n_div)

    def __call__(self):

        from scitbx import matrix
        from random import uniform
        from dials.algorithms.integration import FromBeamVectorToXds
        from dials.algorithms.integration import XdsCoordinateSystem
        from scitbx.array_family import flex
        from profile.tst_profile_helpers import gaussian

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        # Get random x, y, z
        x = uniform(300, 1800)
        y = uniform(300, 1800)
        z = uniform(-10, 0)

        # Get random s1, phi, panel
        s1 = matrix.col(self.detector.get_pixel_lab_coord(
            (x, y))).normalize() * s0_length
        phi = self.scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = self.calculate_bbox(s1, phi, panel)
        x0, x1 = bbox[0], bbox[1]
        y0, y1 = bbox[2], bbox[3]
        z0, z1 = bbox[4], bbox[5]

        # Create pixels in a spot
        size = (z1 - z0, y1 - y0, x1 - x0)
        x0 = (z - z0, y - y0, x - x0)
        sx = (1, 1, 1)
#        sx = tuple([x / 4 for x in x0])
        a = 100
        pixels = gaussian(size, a, x0, sx)
        mask = flex.int(flex.grid(pixels.all()), 1)
        print size, x0
        from matplotlib import pylab
        pylab.imshow(pixels.as_numpy_array()[(z1 - z0) / 2])
        pylab.show()
#        for j in range(z1 - z0):
#            image = pixels.as_numpy_array()[j]
#            pylab.imshow(image, vmin=0, vmax=flex.max(pixels))
#            pylab.show()
#            print image
#        image = pixels[int(x0[0]):int(x0[0])+1,:,:]
#        image.reshape(flex.grid(size[1], size[2]))
#        pylab.imshow(image.as_numpy_array())
#        pylab.show()
#        print 1/0
        print mask.all()

        # Transform the grid
        grid = self.transform(pixels, mask, bbox, s1, phi)

        print flex.min(grid), flex.max(grid)

        gs = grid.all()

        image = grid.as_numpy_array()

        for j in range(11):
            pylab.subplot(3, 4, j + 1)
            print image[j]
            pylab.imshow(image[j], vmin=0, vmax=flex.max(grid))
            pylab.contour(image[j], [self.delta_divergence])
        pylab.show()
#
#        image = grid[int(gs[0]/2):int(gs[0]/2)+1,:,:]
#        print image.as_numpy_array()
#        image.reshape(flex.grid(gs[1], gs[2]))
#        pylab.imshow(image.as_numpy_array())
#        pylab.show()

        print 'OK'


class Test(object):

    def __init__(self):

        import libtbx.load_env
        try:
            dials_regression = libtbx.env.dist_path( 'dials_regression' )
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        import os

        filename = os.path.join(dials_regression,
            'centroid_test_data', 'sweep.json')

        self.tst_reciprocal_space_e3_fraction = \
            TestReciprocalSpaceTransformE3Fraction(filename)
        self.tst_reciprocal_space_transform_detector_lab_coords = \
            TestReciprocalSpaceTransformDetectorLabCoords(filename)
        self.tst_reciprocal_space_transform = \
            TestReciprocalSpaceTransform(filename)

    def run(self):

        #self.tst_reciprocal_space_e3_fraction()
        #self.tst_reciprocal_space_transform_detector_lab_coords()
        self.tst_reciprocal_space_transform()


if __name__ == '__main__':
    test = Test()
    test.run()
