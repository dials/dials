
class TestSpotPredictor:

    def __init__(self):
        from scitbx import matrix
        from scitbx.array_family import flex
        from dials.algorithms import spot_prediction
        from dials_jmp.io import xdsio
        from math import ceil, pi
        from os.path import realpath, dirname, normpath, join
        from dials.model.experiment import Beam, FlatPanelDetector, Goniometer, Scan

        # The XDS files to read from
        test_path = dirname(dirname(dirname(realpath(__file__))))
        integrate_filename = join(test_path, 'data/sim_mx/INTEGRATE.HKL')
        gxparm_filename = join(test_path, 'data/sim_mx/GXPARM.XDS')

        # Read the XDS files
        self.integrate_handle = xdsio.IntegrateFile()
        self.integrate_handle.read_file(integrate_filename)
        self.gxparm_handle = xdsio.GxParmFile()
        self.gxparm_handle.read_file(gxparm_filename)

        # Get the parameters we need from the GXPARM file
        beam = self.gxparm_handle.get_beam()
        gonio = self.gxparm_handle.get_goniometer()
        detector = self.gxparm_handle.get_detector()
        self.ub_matrix = self.gxparm_handle.get_ub_matrix()
        self.ub_matrix = tuple(matrix.sqr(self.ub_matrix).inverse())
        self.unit_cell = self.gxparm_handle.get_unit_cell()
        self.space_group_type = self.gxparm_handle.get_space_group_type()

        # Get the minimum resolution in the integrate file
        d = [self.unit_cell.d(h) for h in self.integrate_handle.hkl]
        self.d_min = min(d)

        # Get the number of frames from the max z value
        xcal, ycal, zcal = zip(*self.integrate_handle.xyzcal)
        gonio.num_frames = int(ceil(max(zcal)))

        detector_d3 = ((matrix.col(detector.x_axis).normalize() *
                        detector.pixel_size[0]) * (0 - detector.origin[0]) +
                       (matrix.col(detector.y_axis).normalize() *
                        detector.pixel_size[1]) * (0 - detector.origin[1]) +
                       detector.distance * matrix.col(detector.normal).normalize())

        self.beam = Beam(beam.direction)
        self.scan = Scan((gonio.starting_frame,
                          gonio.starting_frame + gonio.num_frames),
                         gonio.starting_angle,
                         gonio.oscillation_range, deg=True)
        self.gonio = Goniometer(gonio.rotation_axis)
        self.detector = FlatPanelDetector("",
                            detector.x_axis,
                            detector.y_axis,
                            detector_d3,
                            detector.pixel_size,
                            detector.size,
                            detector.distance,
                            (0, 0))

        # Create the spot predictor
        self.predict_spots = spot_prediction.FlatPanelDetector(
                                           self.beam,
                                           self.detector,
                                           self.gonio,
                                           self.scan,
                                           self.unit_cell,
                                           self.space_group_type,
                                           self.ub_matrix,
                                           self.d_min)

        # Predict the spot locations
        self.reflections = self.predict_spots()

    def test_miller_index_set(self):
        """Ensure we have the whole set of miller indices"""
        gen_hkl = {}
        for r in self.reflections:
            gen_hkl[r.miller_index] = True
        for hkl in self.integrate_handle.hkl:
            assert(gen_hkl[hkl] == True)

    def test_rotation_angles(self):
        """Ensure the rotation angles agree with XDS"""
        from scitbx import matrix

        # Create a dict of lists of xy for each hkl
        gen_phi = {}
        for r in self.reflections:
            hkl = r.miller_index
            phi = r.rotation_angle
            try:
                a = gen_phi[hkl]
                a.append(phi)
                gen_phi[hkl] = a
            except KeyError:
                gen_phi[hkl] = [phi]

        for hkl, xyz in zip(self.integrate_handle.hkl,
                            self.integrate_handle.xyzcal):

            xds_phi = self.scan.starting_angle + xyz[2] * self.scan.oscillation_range

            # Select the nearest xy to use if there are 2
            my_phi = gen_phi[hkl]
            if len(my_phi) == 2:
                my_phi0 = my_phi[0]
                my_phi1 = my_phi[1]
                diff0 = abs(xds_phi - my_phi0)
                diff1 = abs(xds_phi - my_phi1)
                if (diff0 < diff1):
                    my_phi = my_phi0
                else:
                    my_phi = my_phi1
            else:
                my_phi = my_phi[0]

            assert(abs(xds_phi - my_phi) < 0.1)

    def test_beam_vectors(self):
        """Ensure |s1| == |s0|"""
        from scitbx import matrix
        s0_length = matrix.col(self.beam.direction).length()
        for r in self.reflections:
            s1 = r.beam_vector
            s1_length = matrix.col(s1).length()
            assert(abs(s0_length - s1_length) < 1e-7)

    def test_image_coordinates(self):
        """Ensure the image coordinates agree with XDS"""
        from scitbx import matrix

        # Create a dict of lists of xy for each hkl
        gen_xy = {}
        for r in self.reflections:
            hkl = r.miller_index
            xy  = r.image_coord
            try:
                a = gen_xy[hkl]
                a.append(xy)
                gen_xy[hkl] = a
            except KeyError:
                gen_xy[hkl] = [xy]

        for hkl, xyz in zip(self.integrate_handle.hkl,
                            self.integrate_handle.xyzcal):

            xds_xy = (xyz[0], xyz[1])

            # Select the nearest xy to use if there are 2
            my_xy = gen_xy[hkl]
            if len(my_xy) == 2:
                my_xy0 = my_xy[0]
                my_xy1 = my_xy[1]
                diff0 = (matrix.col(xds_xy) - matrix.col(my_xy0)).length()
                diff1 = (matrix.col(xds_xy) - matrix.col(my_xy1)).length()
                if (diff0 < diff1):
                    my_xy = my_xy0
                else:
                    my_xy = my_xy1
            else:
                my_xy = my_xy[0]

            if (abs(xds_xy[0] - my_xy[0]) > 0.1 or
                abs(xds_xy[1] - my_xy[1]) > 0.1):
                print xds_xy, gen_xy[hkl]
            assert(abs(xds_xy[0] - my_xy[0]) < 0.1)
            assert(abs(xds_xy[1] - my_xy[1]) < 0.1)

    def run(self):
        self.test_miller_index_set()
        self.test_rotation_angles()
        self.test_beam_vectors()
        self.test_image_coordinates()
        print "OK"

if __name__ == '__main__':
    test = TestSpotPredictor()
    test.run()
