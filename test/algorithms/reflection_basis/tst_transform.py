
class TestForward(object):
    def __init__(self, filename):
        from dials.model.serialize import load
        from dials.algorithms.reflection_basis import transform
        from dials.algorithms.integration import BBoxCalculator
        from math import pi

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        print self.detector.get_d_matrix()
        print self.gonio.get_rotation_axis()
        print self.beam.get_direction()
        print self.detector.get_beam_centre(self.beam.get_s0())
        self.beam.set_direction((0.0, 0.0, 1.0))
        self.gonio.set_rotation_axis((1.0, 0.0, 0.0))
        self.detector.set_frame((1.0, 0.0, 0.0),
                                (0.0, 1.0, 0.0),
                                (-150, -150, -200))

        # Set some parameters
        self.sigma_divergence =self.beam.get_sigma_divergence(deg=False)
        self.mosaicity = 0.157 * pi / 180
        self.n_sigma = 3
        self.grid_size = 7
        self.delta_divergence = self.n_sigma * self.sigma_divergence

        step_size = self.delta_divergence / self.grid_size
        self.delta_divergence2 = self.delta_divergence + step_size * 0.5
        self.delta_mosaicity = self.n_sigma * self.mosaicity

        # Create the bounding box calculator
        self.calculate_bbox = BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            self.delta_divergence2,
            self.delta_mosaicity)

        # Initialise the transform
        self.transform = transform.Forward(
            self.beam, self.detector, self.scan,
            self.mosaicity, self.n_sigma, self.grid_size)

    def __call__(self):

        from scitbx import matrix
        from random import uniform
        from dials.algorithms.reflection_basis import CoordinateSystem
        from scitbx.array_family import flex
        from time import time

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        # Get random x, y, z
        x = uniform(300, 1800)
        y = uniform(300, 1800)
        z = uniform(-10, 0)
        print x, y, z

#        x, y, z = 500, 500, 10
        x, y, z = 993.376692565, 1272.1175661, -7.9419902089

        # Get random s1, phi, panel
        s1 = matrix.col(self.detector.get_pixel_lab_coord(
            (x, y))).normalize() * s0_length
        phi = self.scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = self.calculate_bbox(s1, phi, panel)
        x0, x1, y0, y1, z0, z1 = bbox

        from dials.algorithms.reflection_basis import FromDetector


        print bbox, x, y, z

        # Create the coordinate system
        cs = CoordinateSystem(m2, s0, s1, phi)

        from_detector = FromDetector(cs, self.detector.get_d_matrix())
        step_size = self.delta_divergence / self.grid_size
        e0, e1 = from_detector(self.detector.pixel_to_millimeter((x, y)))
        e00, e01 = from_detector(self.detector.pixel_to_millimeter((x0, y0)))
        e10, e11 = from_detector(self.detector.pixel_to_millimeter((x0, y1)))
        e20, e21 = from_detector(self.detector.pixel_to_millimeter((x1, y0)))
        e30, e31 = from_detector(self.detector.pixel_to_millimeter((x1, y1)))
        offset = self.grid_size + 0.5
        print e0 / step_size + offset, e1 / step_size + offset
        print e00 / step_size + offset, e01 / step_size + offset
        print e10 / step_size + offset, e11 / step_size + offset
        print e20 / step_size + offset, e21 / step_size + offset
        print e30 / step_size + offset, e31 / step_size + offset
        print self.delta_divergence, self.delta_divergence2,

        # Create the image
        image = flex.double(flex.grid(z1 - z0, y1 - y0, x1 - x0), 1)
        mask = flex.bool(flex.grid(image.all()), True)

        # Transform the image to the grid
        grid = self.transform(cs, bbox, image, mask)

        #print flex.min(grid), flex.max(grid)
        print grid.as_numpy_array()[4,:,:]
        print grid.all()

        from dials.algorithms.reflection_basis import transform
        step_size = self.delta_divergence / self.grid_size
        s1_map = transform.beam_vector_map(self.detector, self.beam, True)
        index = transform.GridIndexGenerator(cs, x0, y0, (step_size, step_size),
            self.grid_size, s1_map)

        indices = []
        labels = []
        for j in range(y1 - y0 + 1):
            for i in range(x1 - x0 + 1):
                indices.append(index(j,i))
                labels.append('{0}, {1}'.format(i, j))
        x, y = zip(*indices)

        #x = [xx - 0.5 for xx in x]
        #y = [yy - 0.5 for yy in y]
        #print x, y
        from dials.algorithms.polygon import simple_area
        indices = []
        area = []
        from_detector = FromDetector(cs, self.detector.get_d_matrix())
        #index2 = lambda j, i: matrix.col(from_detector(self.detector.pixel_to_millimeter((i, j)))) / step_size + matrix.col((offset, offset))


        def index(j, i):
            mm = self.detector.pixel_to_millimeter((i + x0, j + y0))
            ee = matrix.col(from_detector(mm))
            return (ee[0] / step_size + offset, ee[1] / step_size + offset)

        for j in range(y1 - y0):
            for i in range(x1 - x0):
                gx00, gy00 = index(j, i)
                gx01, gy01 = index(j, i+1)
                gx10, gy10 = index(j+1,i)
                gx11, gy11 = index(j+1,i+1)
                poly = flex.vec2_double([
                    (gx00, gy00),
                    (gx01, gy01),
                    (gx11, gy11),
                    (gx10, gy10)])
                indices.append((gx00, gy00))
                area.append(simple_area(poly))

        from matplotlib import pylab
        x, y = zip(*indices)
        pylab.plot(y, area)
        pylab.show()
#        for i in range(9):
#            pylab.subplot(3, 3, i)
#            pylab.imshow(grid.as_numpy_array()[i,:,:], origin='bottom')
#            pylab.scatter(x, y)

        #pylab.imshow(grid.as_numpy_array()[4,:,:], origin='bottom', interpolation='none')
#        pylab.scatter(x, y)
#        for j in range(9):
#            for i in range(9):
#                ii = i - 0.5
#                jj = j - 0.5
#                px = [ii, ii + 1, ii + 1, ii, ii]
#                py = [jj, jj, jj + 1, jj + 1, jj]
#                pylab.plot(px, py, color='white')
#        for ll, xx, yy in zip(labels, x, y):
#            pylab.annotate(ll, xy = (xx, yy))
        #pylab.show()

        # Test passed
        print 'OK'


class TestReverse(object):
    def __init__(self, filename):
        from dials.model.serialize import load
        from dials.algorithms.reflection_basis import transform
        from math import pi

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.scan = self.sweep.get_scan()

        # Set some parameters
        self.mosaicity = 0.157 * pi / 180
        self.n_sigma = 3
        self.grid_size = 4

        # Initialise the transform
        self.transform = transform.Reverse(
            self.beam, self.detector, self.scan,
            self.mosaicity, self.n_sigma, self.grid_size)

    def __call__(self):

        # Test passed
        print 'OK'

class Test(object):
    def __init__(self):

        import os
        import libtbx.load_env

        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        # Set the sweep filename and load the sweep
        filename = os.path.join(
            dials_regression,
            'centroid_test_data',
            'sweep.json')

        self.tst_forward = TestForward(filename)
        self.tst_reverse = TestReverse(filename)

    def run(self):
        self.tst_forward()
        #self.tst_reverse()


if __name__ == '__main__':
    test = Test()
    test.run()
