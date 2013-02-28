from __future__ import division

class TestRayPredictor:

    def __init__(self):
        from scitbx import matrix
        from dials.algorithms.spot_prediction import RayPredictor
        from dials.algorithms.spot_prediction import IndexGenerator
        from iotbx.xds import xparm, integrate_hkl
        from dials.util import io
        from math import ceil
        from os.path import realpath, dirname, join

        # The XDS files to read from
        test_path = dirname(dirname(dirname(realpath(__file__))))
        integrate_filename = join(test_path, 'data/sim_mx/INTEGRATE.HKL')
        gxparm_filename = join(test_path, 'data/sim_mx/GXPARM.XDS')

        # Read the XDS files
        self.integrate_handle = integrate_hkl.reader()
        self.integrate_handle.read_file(integrate_filename)
        self.gxparm_handle = xparm.reader()
        self.gxparm_handle.read_file(gxparm_filename)        

        # Get the parameters we need from the GXPARM file
        models = io.read_models_from_file(gxparm_filename)
        self.beam = models.get_beam()
        self.gonio = models.get_goniometer()
        self.detector = models.get_detector()
        self.scan = models.get_scan()        

        # Get crystal parameters
        self.ub_matrix = io.get_ub_matrix_from_xparm(self.gxparm_handle)
        self.unit_cell = io.get_unit_cell_from_xparm(self.gxparm_handle)
        self.space_group_type = io.get_space_group_type_from_xparm(self.gxparm_handle)

        # Get the minimum resolution in the integrate file
        d = [self.unit_cell.d(h) for h in self.integrate_handle.hkl]
        self.d_min = min(d)

        # Get the number of frames from the max z value
        xcal, ycal, zcal = zip(*self.integrate_handle.xyzcal)
        self.scan.image_range = (self.scan.image_range[0],
                                self.scan.image_range[0] + int(ceil(max(zcal))))

        # Print stuff
#        print self.beam
#        print self.gonio
#        print self.detector
#        print self.scan

        # Create the index generator
        self.generate_indices = IndexGenerator(self.unit_cell,
            self.space_group_type, self.d_min)

        s0 = self.beam.direction
        m2 = self.gonio.rotation_axis
        UB = self.ub_matrix
        dphi = self.scan.get_oscillation_range()

        # Create the ray predictor
        self.predict_rays = RayPredictor(s0, m2, UB, dphi)

        # Predict the spot locations
        self.reflections = self.predict_rays(self.generate_indices.to_array())

    def test_miller_index_set(self):
        """Ensure we have the whole set of miller indices"""
        gen_hkl = {}
        print len(self.reflections)
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

        # For each hkl in the xds file
        for hkl, xyz in zip(self.integrate_handle.hkl,
                            self.integrate_handle.xyzcal):

            # Calculate the XDS phi value
            xds_phi = self.scan.oscillation[0] + xyz[2]*self.scan.oscillation[1]

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

            # Check the Phi values are the same
            assert(abs(xds_phi - my_phi) < 0.1)

        print "OK"

    def test_beam_vectors(self):
        """Ensure |s1| == |s0|"""
        from scitbx import matrix
        s0_length = matrix.col(self.beam.direction).length()
        for r in self.reflections:
            s1 = r.beam_vector
            s1_length = matrix.col(s1).length()
            assert(abs(s0_length - s1_length) < 1e-7)

        print "OK"

    def run(self):
        self.test_miller_index_set()
        self.test_rotation_angles()
        self.test_beam_vectors()

if __name__ == '__main__':
    test = TestRayPredictor()
    test.run()
