
class Test(object):

    def __init__(self):
        from dials.algorithms.integration import BBoxCalculator
        from dials.model.serialize import load
        from math import pi

        import libtbx.load_env
        try:
            dials_regression = libtbx.env.dist_path( 'dials_regression' )
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        import os

        filename = os.path.join(dials_regression,
            'centroid_test_data', 'sweep.json')

        sweep = load.sweep(filename)

        # Get the models
        self.beam = sweep.get_beam()
        self.detector = sweep.get_detector()
        self.gonio = sweep.get_goniometer()
        self.scan = sweep.get_scan()

        # Set the delta_divergence/mosaicity
        self.n_sigma = 5
        self.sigma_divergence = 0.060 * pi / 180
        self.mosaicity = 0.154 * pi / 180
        self.delta_divergence = self.n_sigma * self.sigma_divergence
        self.delta_mosaicity = self.n_sigma * self.mosaicity

        # Create the bounding box calculator
        self.calculate_bbox = BBoxCalculator(
            self.beam, self.detector, self.gonio, self.scan,
            self.delta_divergence,
            self.delta_mosaicity)

    def run(self):
        self.tst_outer_bounds()
        self.tst_radius()

    def tst_outer_bounds(self):

        from scitbx import matrix
        from random import uniform
        from dials.algorithms import reflection_basis as rb

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        for i in range(1000):

            # Get random x, y, z
            x = uniform(0, 2000)
            y = uniform(0, 2000)
            z = uniform(0, 1000)

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

            # Calculate the rotation angle for each point
            phi_dash1 = self.scan.get_angle_from_array_index(z1, deg = False)
            phi_dash2 = self.scan.get_angle_from_array_index(z2, deg = False)

            # Create the XDS coordinate system
            xcs = rb.CoordinateSystem(m2, s0, s1, phi)

            # Calculate reciprocal space coordinates at each point
            to_reciprocal_space = rb.FromBeamVectorAndRotationAngle(xcs)
            e11, e21, e31 = to_reciprocal_space(s1, phi_dash1)
            e12, e22, e32 = to_reciprocal_space(s1, phi_dash2)

            # Check vertical edges
            for j in range(bbox[2], bbox[3] + 1):
                xyz1 = self.detector.get_pixel_lab_coord((bbox[0], j))
                xyz2 = self.detector.get_pixel_lab_coord((bbox[1] + 1, j))
                sdash1 = matrix.col(xyz1).normalize() * s0_length
                sdash2 = matrix.col(xyz2).normalize() * s0_length
                e11, e21, e3 = to_reciprocal_space(sdash1, phi)
                e12, e22, e3 = to_reciprocal_space(sdash2, phi)
                assert(abs(e11) >= self.delta_divergence or
                       abs(e21) >= self.delta_divergence)
                assert(abs(e12) >= self.delta_divergence or
                       abs(e22) >= self.delta_divergence)

            # Check horizontal edges
            for i in range(bbox[0], bbox[1] + 1):
                xyz1 = self.detector.get_pixel_lab_coord((i, bbox[2]))
                xyz2 = self.detector.get_pixel_lab_coord((i, bbox[3] + 1))
                sdash1 = matrix.col(xyz1).normalize() * s0_length
                sdash2 = matrix.col(xyz2).normalize() * s0_length
                e11, e21, e3 = to_reciprocal_space(sdash1, phi)
                e12, e22, e3 = to_reciprocal_space(sdash2, phi)
                assert(abs(e11) >= self.delta_divergence or
                       abs(e21) >= self.delta_divergence)
                assert(abs(e12) >= self.delta_divergence or
                       abs(e22) >= self.delta_divergence)

            # All e3 coords >= delta_mosaicity
            assert(abs(e31) >= self.delta_mosaicity)
            assert(abs(e32) >= self.delta_mosaicity)

        print 'OK'

    def tst_radius(self):
        from scitbx import matrix
        from random import uniform
        from dials.algorithms import reflection_basis as rb
        from math import sqrt

        s0 = self.beam.get_s0()
        m2 = self.gonio.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()

        radius12 = self.delta_divergence
        radius3 = self.delta_mosaicity

        for i in range(1000):

            # Get random x, y, z
            x = uniform(0, 2000)
            y = uniform(0, 2000)
            z = uniform(0, 1000)

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

            # Calculate the rotation angle for each point
            phi_dash1 = self.scan.get_angle_from_array_index(z1, deg = False)
            phi_dash2 = self.scan.get_angle_from_array_index(z2, deg = False)

            # Create the XDS coordinate system
            xcs = rb.CoordinateSystem(m2, s0, s1, phi)

            # Calculate reciprocal space coordinates at each point
            to_reciprocal_space = rb.FromBeamVectorAndRotationAngle(xcs)

            # Check vertical edges
            for j in range(bbox[2], bbox[3] + 1):
                xyz1 = self.detector.get_pixel_lab_coord((bbox[0], j))
                xyz2 = self.detector.get_pixel_lab_coord((bbox[1] + 1, j))
                sdash1 = matrix.col(xyz1).normalize() * s0_length
                sdash2 = matrix.col(xyz2).normalize() * s0_length
                e11, e21, e31 = to_reciprocal_space(sdash1, phi)
                e12, e22, e31 = to_reciprocal_space(sdash2, phi)
                assert(sqrt(e11**2 + e21**2) >= radius12)
                assert(sqrt(e12**2 + e22**2) >= radius12)

            # Check horizontal edges
            for i in range(bbox[0], bbox[1] + 1):
                xyz1 = self.detector.get_pixel_lab_coord((i, bbox[2]))
                xyz2 = self.detector.get_pixel_lab_coord((i, bbox[3] + 1))
                sdash1 = matrix.col(xyz1).normalize() * s0_length
                sdash2 = matrix.col(xyz2).normalize() * s0_length
                e11, e21, e32 = to_reciprocal_space(sdash1, phi)
                e12, e22, e32 = to_reciprocal_space(sdash2, phi)
                assert(sqrt(e11**2 + e21**2) >= radius12)
                assert(sqrt(e12**2 + e22**2) >= radius12)

        print 'OK'

if __name__ == '__main__':
    test = Test()
    test.run()
