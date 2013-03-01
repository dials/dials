from __future__ import division

def run():
    from scitbx import matrix
    from iotbx.xds import xparm, integrate_hkl
    from dials.util import io
    from math import ceil
    from dials.algorithms.spot_prediction import RotationAngles
    from os.path import realpath, dirname, join
    import dxtbx

    # The XDS files to read from
    test_path = dirname(dirname(dirname(realpath(__file__))))
    integrate_filename = join(test_path, 'data/sim_mx/INTEGRATE.HKL')
    gxparm_filename = join(test_path, 'data/sim_mx/GXPARM.XDS')

    # Read the XDS files
    integrate_handle = integrate_hkl.reader()
    integrate_handle.read_file(integrate_filename)
    gxparm_handle = xparm.reader()
    gxparm_handle.read_file(gxparm_filename)

    # Get the parameters we need from the GXPARM file
    models = dxtbx.load(gxparm_filename)
    beam = models.get_beam()
    gonio = models.get_goniometer()
    detector = models.get_detector()
    scan = models.get_scan()
    
    # Get the crystal parameters
    ub_matrix = io.get_ub_matrix_from_xparm(gxparm_handle)
    unit_cell = io.get_unit_cell_from_xparm(gxparm_handle)
    space_group_type = io.get_space_group_type_from_xparm(gxparm_handle)

    # Get the minimum resolution in the integrate file
    d = [unit_cell.d(h) for h in integrate_handle.hkl]
    d_min = min(d)

    # Get the number of frames from the max z value
    xcal, ycal, zcal = zip(*integrate_handle.xyzcal)
    num_frames = int(ceil(max(zcal)))
    scan.image_range = (scan.image_range[0], 
                        scan.image_range[0] + num_frames - 1)
    
    # Create the rotation angle object
    ra = RotationAngles(beam.direction, gonio.rotation_axis)

    # Setup the matrices
    ub = matrix.sqr(ub_matrix)
    s0 = matrix.col(beam.direction)
    m2 = matrix.col(gonio.rotation_axis)

    # For all the miller indices
    for h in integrate_handle.hkl:
      h = matrix.col(h)

      # Calculate the angles
      angles = ra(h, ub)

      # For all the angles
      for phi in angles:
          r = m2.axis_and_angle_as_r3_rotation_matrix(angle=phi)
          pstar = r * ub * h
          s1 = s0 + pstar
          assert(abs(s1.length() - s0.length()) < 1e-7)

    # Test Passed
    print "OK"

if __name__ == '__main__':
    run()
