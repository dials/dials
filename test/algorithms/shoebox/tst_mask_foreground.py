from __future__ import division

class Test(object):

    def __init__(self):
        import os
        import libtbx.load_env
        from dials.model.serialize import load
        from dials.algorithms.shoebox import BBoxCalculator
        from dials.algorithms.shoebox import MaskForeground

        try:
            dials_regression = libtbx.env.dist_path('dials_regression')
        except KeyError, e:
            print 'FAIL: dials_regression not configured'
            return

        # Set the sweep filename and load the sweep
        sweep_filename = os.path.join(dials_regression, 'centroid_test_data',
            'sweep.json')
        crystal_filename = os.path.join(dials_regression, 'centroid_test_data',
            'crystal.json')

        # Load the sweep
        self.sweep = load.sweep(sweep_filename)
        self.crystal = load.crystal(crystal_filename)
        self.beam = self.sweep.get_beam()
        self.detector = self.sweep.get_detector()
        self.goniometer = self.sweep.get_goniometer()
        self.scan = self.sweep.get_scan()
        self.delta_d = 5 * self.beam.get_sigma_divergence(deg=False)
        self.delta_m = 5 * self.crystal.get_mosaicity(deg=False)

        # Get the function object to mask the foreground
        self.mask_foreground = MaskForeground(self.beam, self.detector,
            self.goniometer, self.scan, self.delta_d, self.delta_m)

        # Get the function object to calcyulate the bounding box
        self.calculate_bbox = BBoxCalculator(self.beam, self.detector,
            self.goniometer, self.scan, self.delta_d, self.delta_m)

    def run(self):
        from scitbx import matrix
        from dials.algorithms import shoebox
        from dials.algorithms.reflection_basis import CoordinateSystem

        s0 = self.beam.get_s0()
        m2 = self.goniometer.get_rotation_axis()
        s0_length = matrix.col(self.beam.get_s0()).length()
        width, height = self.detector.get_image_size()

        # Generate some reflections
        reflections = self.generate_reflections(100)

        # Mask the foreground in each
        self.mask_foreground(reflections)

        # Loop through all the reflections and check the mask values
        for r in reflections:
            mask = r.shoebox_mask[0:1,:,:]
            x0, x1, y0, y1, z0, z1 = r.bounding_box
            s1 = r.beam_vector
            phi = r.rotation_angle
            cs = CoordinateSystem(m2, s0, s1, phi)
            for j in range(y1 - y0):
                for i in range(x1 - x0):
                    value1 = mask[0, j, i]
                    s1 = self.detector.get_pixel_lab_coord(
                        (x0 + i + 0.5, y0 + j + 0.5))
                    s1 = matrix.col(s1).normalize() * s0_length
                    e1, e2 = cs.from_beam_vector(s1)
                    aa = (e1 / self.delta_d)**2
                    bb = (e2 / self.delta_d)**2
                    if (x0 + i < 0 or y0 + j < 0 or
                        x0 + i > width or y0 + j > height):
                        value2 = 0
                    else:
                        if (aa + bb <= 1.0):
                            value2 = shoebox.MaskCode.Foreground
                        else:
                            value2 = shoebox.MaskCode.Background
                    assert(value1 == value2)

        # Test passed
        print 'OK'

    def generate_reflections(self, num):
        from dials.model.data import ReflectionList
        from random import randint
        from scitbx import matrix
        from scitbx.array_family import flex
        rlist = ReflectionList(num)
        s0_length = matrix.col(self.beam.get_s0()).length()
        for i in range(num):
            x = randint(0, 2000)
            y = randint(0, 2000)
            z = randint(0, 100)
            s1 = self.detector.get_pixel_lab_coord((x, y))
            s1 = matrix.col(s1).normalize() * s0_length
            phi = self.scan.get_angle_from_array_index(z)
            rlist[i].beam_vector = s1
            rlist[i].rotation_angle = phi
        self.calculate_bbox(rlist)
        for i in range(num):
            bbox = rlist[i].bounding_box
            xsize = bbox[1] - bbox[0]
            ysize = bbox[3] - bbox[2]
            zsize = bbox[5] - bbox[4]
            rlist[i].shoebox_mask = flex.int(flex.grid(zsize, ysize, xsize))
        return rlist

if __name__ == '__main__':
    test = Test()
    test.run()
