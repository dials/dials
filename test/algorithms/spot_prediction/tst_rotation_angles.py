
def run():
    from scitbx import matrix
    from scitbx.array_family import flex
    from dials_jmp.io import xdsio
    from math import ceil, pi
    from dials.algorithms.spot_prediction import RotationAngles
    from os.path import realpath, dirname, normpath, join

    # The XDS files to read from
    test_path = dirname(dirname(dirname(realpath(__file__))))
    integrate_filename = join(test_path, 'data/sim_mx/INTEGRATE.HKL')
    gxparm_filename = join(test_path, 'data/sim_mx/GXPARM.XDS')

    # Read the XDS files
    integrate_handle = xdsio.IntegrateFile()
    integrate_handle.read_file(integrate_filename)
    gxparm_handle = xdsio.GxParmFile()
    gxparm_handle.read_file(gxparm_filename)

    # Get the parameters we need from the GXPARM file
    beam = gxparm_handle.get_beam()
    gonio = gxparm_handle.get_goniometer()
    detector = gxparm_handle.get_detector()
    ub_matrix = gxparm_handle.get_ub_matrix()
    ub_matrix = tuple(matrix.sqr(ub_matrix).inverse())
    unit_cell = gxparm_handle.get_unit_cell()
    space_group_type = gxparm_handle.get_space_group_type()

    # Get the minimum resolution in the integrate file
    d = [unit_cell.d(h) for h in integrate_handle.hkl]
    d_min = min(d)

    # Get the number of frames from the max z value
    xcal, ycal, zcal = zip(*integrate_handle.xyzcal)
    gonio.num_frames = int(ceil(max(zcal)))

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
