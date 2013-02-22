import unittest

class TestIndexGenerator(unittest.TestCase):

    def setUp(self):
        pass

    def test_range(self):

        from dials.io import xdsio
        from dials.spot_prediction import IndexGenerator
        import numpy

        # The XDS files to read from
        integrate_filename = './test/data/sim_mx/INTEGRATE.HKL'
        gxparm_filename = './test/data/sim_mx/GXPARM.XDS'

        # Read the XDS files
        integrate_handle = xdsio.IntegrateFile()
        integrate_handle.read_file(integrate_filename)
        gxparm_handle = xdsio.GxParmFile()
        gxparm_handle.read_file(gxparm_filename)

        # Get the parameters we need from the GXPARM file
        d_min = 1.6
        unit_cell = gxparm_handle.get_unit_cell()
        space_group_type = gxparm_handle.get_space_group_type()

        # Generate the indices
        index_generator = IndexGenerator(unit_cell, space_group_type, True, d_min)
        miller_indices = index_generator.to_array()

        # Get individual generated hkl
        gen_h = [hkl[0] for hkl in miller_indices]
        gen_k = [hkl[1] for hkl in miller_indices]
        gen_l = [hkl[2] for hkl in miller_indices]

        # Get individual xds generated hkl
        xds_h = [hkl[0] for hkl in integrate_handle.hkl]
        xds_k = [hkl[0] for hkl in integrate_handle.hkl]
        xds_l = [hkl[0] for hkl in integrate_handle.hkl]

        # Get min/max generated hkl
        min_gen_h, max_gen_h = numpy.min(gen_h), numpy.max(gen_h)
        min_gen_k, max_gen_k = numpy.min(gen_k), numpy.max(gen_k)
        min_gen_l, max_gen_l = numpy.min(gen_l), numpy.max(gen_l)

        # Get min/max xds generated hkl
        min_xds_h, max_xds_h = numpy.min(xds_h), numpy.max(xds_h)
        min_xds_k, max_xds_k = numpy.min(xds_k), numpy.max(xds_k)
        min_xds_l, max_xds_l = numpy.min(xds_l), numpy.max(xds_l)

        # Ensure we have the whole xds range  in the generated set
        self.assertTrue(min_gen_h <= min_xds_h and max_gen_h >= max_xds_h)
        self.assertTrue(min_gen_k <= min_xds_k and max_gen_k >= max_xds_k)
        self.assertTrue(min_gen_l <= min_xds_l and max_gen_l >= max_xds_l)

        # Ensure we have the whole xds range (takes ages)
#        for hkl in integrate_handle.hkl:
#            self.assertTrue(hkl in miller_indices)


if __name__ == '__main__':
    unittest.main()
