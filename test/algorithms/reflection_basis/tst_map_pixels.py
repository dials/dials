from __future__ import division

def evaluate_gaussian(x, a, x0, sx):

    from math import exp

    assert(len(x) == len(x0))
    assert(len(x) == len(sx))

    g = 0.0
    for xi, x0i, sxi in zip(x, x0, sx):
        g += (xi - x0i)**2 / (2.0 * sxi**2)
    return a * exp(-g)

def gaussian(size, a, x0, sx):

    from scitbx.array_family import flex

    result = flex.double(flex.grid(size))
    index = [0 for i in range(len(size))]
    while True:
        result[index[::-1]] = evaluate_gaussian(index[::-1], a, x0, sx)
        for j in range(len(size)):
            index[j] += 1
            if index[j] < size[::-1][j]:
                break
            index[j] = 0
            if j == len(size) - 1:
                return result

class TestForward(object):
    def __init__(self, filename):
        from dials.model.serialize import load
        from dials.algorithms.reflection_basis import transform
        from dials.algorithms.shoebox import BBoxCalculator
        from math import pi

        # Load the sweep
        self.sweep = load.sweep(filename)

        # Get the models
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.gonio = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()

        # Set some parameters
        self.sigma_divergence =self.beam.get_sigma_divergence(deg=False)
        self.mosaicity = 0.157 * pi / 180
        self.n_sigma = 3
        self.grid_size = 4
        self.delta_divergence = self.n_sigma * self.sigma_divergence

        step_size = self.delta_divergence / self.grid_size
        self.delta_divergence2 = self.delta_divergence + step_size * 0.5
        self.delta_mosaicity = self.n_sigma * self.mosaicity

        # Create the bounding box calculator
        self.calculate_bbox = BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            self.delta_divergence2,
            self.delta_mosaicity)

        self.s1_map = transform.beam_vector_map(self.detector, self.beam, True)

        # Initialise the transform
        self.map_pixels = transform.MapPixelsForward(
            self.s1_map, self.grid_size, (step_size, step_size))

        self.map_frames = transform.MapFramesForward(
            self.scan.get_oscillation()[0],
            self.scan.get_oscillation()[1],
            self.mosaicity,
            self.n_sigma,
            self.grid_size)

    def __call__(self):
        self.tst_output()

    def tst_output(self):

        from dials.algorithms.reflection_basis import CoordinateSystem
        from dials.algorithms.reflection_basis import transform
        from scitbx import matrix
        from scitbx.array_family import flex

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        # Get random x, y, z
        x = 1000
        y = 1000
        z = 0

        # Get random s1, phi, panel
        s1 = matrix.col(self.detector.get_pixel_lab_coord(
            (x, y))).normalize() * s0_length
        phi = self.scan.get_angle_from_array_index(z, deg=False)
        panel = 0

        # Calculate the bounding box
        bbox = self.calculate_bbox(s1, phi, panel)
        x0, x1, y0, y1, z0, z1 = bbox

        # Create the coordinate system
        cs = CoordinateSystem(m2, s0, s1, phi)

        # Create the image
        #image = flex.double(flex.grid(z1 - z0, y1 - y0, x1 - x0), 1)
        image = gaussian((z1 - z0, y1 - y0, x1 - x0), 10.0,
            (z - z0, y - y0, x - x0), (2.0, 2.0, 2.0))
        mask = flex.bool(flex.grid(image.all()), True)

        z_fraction = self.map_frames(bbox[4:], phi, cs.zeta())

        # Transform the image to the grid
        output = self.map_pixels(cs, bbox, image, mask, z_fraction)
        import pickle
        import os
        expected = pickle.load(open(os.path.join(
            os.path.dirname(__file__), 'expected.p'), 'r'))

        eps = 1e-7
        assert(len(expected) == len(output))
        for e, o in zip(expected, output):
            assert(abs(e - o) < eps)

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

    def run(self):
        self.tst_forward()


if __name__ == '__main__':
    test = Test()
    test.run()
