from scitbx import matrix


class GxParmFile:
    """A class to read the GXPARM.XDS file used in XDS"""
    
    def __init__(self):
        pass

    def read_file(self, filename):
        """Read the GXPARAM.XDS file.
        
        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.
        
        :param filename: The path to the file

        """
        # Read the text from the file and split into an array of tokens
        lines = open(filename, 'r').readlines()
        tokens = [map(float, l.split()) for l in lines]
        
        # Read the parameters from the list of tokens
        self.starting_frame    = tokens[0][0]
        self.starting_angle    = tokens[0][1]
        self.oscillation_range = tokens[0][2]
        self.rotation_axis     = tuple(tokens[0][3:6])
        self.wavelength        = tokens[1][0]
        self.beam_vector       = tuple(tokens[1][1:4])
        self.detector_size     = tuple(map(int, tokens[2][0:2]))
        self.pixel_size        = tuple(tokens[2][2:4])
        self.detector_distance = tokens[3][0]
        self.detector_origin   = tuple(tokens[3][1:3])
        self.detector_x_axis   = tuple(tokens[4])
        self.detector_y_axis   = tuple(tokens[5])
        self.detector_normal   = tuple(tokens[6])
        self.space_group       = int(tokens[7][0])
        self.unit_cell         = tuple(tokens[7][1:7])
        self.unit_cell_a_axis  = tuple(tokens[8])
        self.unit_cell_b_axis  = tuple(tokens[9])
        self.unit_cell_c_axis  = tuple(tokens[10])
        
        
class IntegrateFile:
    """A class to read the INTEGRATE.HKL file used in XDS"""

    def __init__(self):
        """Initialise the file contents."""
    
        self._header = {}
    
        self.space_group = None
        self.unit_cell = None
        self.detector_size = None
        self.pixel_size = None
        self.starting_frame = None
        self.starting_angle = None
        self.oscillation_range = None
        self.rotation_axis = None
        self.wavelength = None
        self.beam_vector = None
        self.detector_x_axis = None
        self.detector_y_axis = None
        self.detector_origin = None
        self.detector_origin = None
        self.detector_distance = None
        self.unit_cell_a_axis = None
        self.unit_cell_b_axis = None
        self.unit_cell_c_axis = None
    
        self.hkl = []
        self.iobs = []
        self.sigma = []
        self.xyzcal = []
        self.rlp = []
        self.peak = []
        self.corr = []
        self.maxc = []
        self.xyzobs = []
        self.alfbet0 = []
        self.alfbet1 = []
        self.psi = []
    
    
    def read_file(self, filename):
        """Read the INTEGRATE.HKL file.
        
        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.
        
        :param filename: The path to the file
        
        """
        # Read the lines from the file
        lines = open(filename, 'r').readlines()

        # Loop through the lines in the file. First off, parse the header
        # lines until we reach !END_OF_HEADER. Then parse the data lines
        # until we read !END_OF_DATA
        in_header = True
        for l in lines:
            if in_header:
                if l.strip().startswith('!END_OF_HEADER'):
                    in_header = False
                    continue
                else:
                    pass
                    #self._parse_header_line(l)
            else:
                if l.strip().startswith('!END_OF_DATA'):
                    break
                else:
                    self._parse_data_line(l)

      
    def _parse_data_line(self, line):
        """Parse a data line from the Integrate.hkl file
        
        :param line: The line to parse
        
        """
        # Split the tokens
        tokens = line.split()
        tokens = map(int, tokens[0:3]) + map(float, tokens[3:])

        # Get the reflection information and append to the lists
        self.hkl.append    (tuple(tokens[0:3]))
        self.iobs.append   (tokens[3])
        self.sigma.append  (tokens[4])
        self.xyzcal.append (tuple(tokens[5:8]))
        self.rlp.append    (tokens[8])
        self.peak.append   (tokens[9])
        self.corr.append   (tokens[10])
        self.maxc.append   (tokens[11])
        self.xyzobs.append (tuple(tokens[12:15]))
        self.alfbet0.append(tuple(tokens[15:17]))
        self.alfbet1.append(tuple(tokens[17:19]))
        self.psi.append    (tokens[19])
  
  
def calculate_reflection_detector_coordinates(indices, s0, m2, b1, b2, b3, x0, 
                                              y0, f, d1, d2, d3, phi_list):
    """Calculate the pixel coordinates of each reflection.
    
    Calculate the diffracted beam wave vector and find where it intersects
    with the dectector plane.
    
    :param indices: The miller indices
    :param s0: The incident beam wave vector
    :param m2: The crystal rotation axis
    :param b1: The unit cell a axis unit vector
    :param b2: The unit cell b axis unit vector
    :param b3: The unit cell c axis unit vector
    :param x0: The x coordinate of the detector origin
    :param y0: The y coordinate of the detector origin
    :param f: The distance from the crystal to the detector
    :param d1: The detector x axis vector
    :param d2: The detector y axis vector
    :param d3: The detector normal unit vector
    :param phi: The rotation angle
    :returns: A list of pixel indices
    
    """
    # Ensure everything is a matrix
    m2 = matrix.col(m2)
    s0 = matrix.col(s0)
    b1 = matrix.col(b1)
    b2 = matrix.col(b2)
    b3 = matrix.col(b3)
    d1 = matrix.col(d1)
    d2 = matrix.col(d2)
    d3 = matrix.col(d3)

    # Calculate the reciprocal basis of the crystal lattice coordinate system
    b1_star = b2.cross(b3) / b1.dot(b2.cross(b3))
    b2_star = b3.cross(b1) / b1.dot(b2.cross(b3))
    b3_star = b1.cross(b2) / b1.dot(b2.cross(b3))
    
    # Loop through all the miller indices. Calculate the pixel coordinates of
    # the reflection and append them to the coordinate list.
    coords = []
    for phi, (h, k, l) in zip(phi_list, indices):
 
        # Calculate the reciprocal lattice vector and rotate it by phi
        # about the rotation axis. Then calculate the diffracted beam vector, s.
        p_star0 = h*b1_star + k*b2_star + l*b3_star
        p_star = p_star0.rotate(axis=m2, angle=phi, deg=True)
        s = s0 + p_star

        # Calculate the point at which the diffraction beam vector intersects 
        # the detector plane. (only when s.d3 > 0)
        s_dot_d3 = s.dot(d3)
        if (s_dot_d3 > 0):
            x = x0 + f * s.dot(d1) / s_dot_d3
            y = y0 + f * s.dot(d2) / s_dot_d3
            coords.append((x, y))

    # Return the pixel coordinates
    return coords


def validate_calculated_coordinates(xcal, ycal, xexp, yexp):
    """Check the calculated coordinates against the expected"""

    from math import sqrt
    import numpy

    # Print a message if the wrong length
    if len(xcal) != len(xexp):
        print "len(coords) doesn't match len(xcal)"

    # Calculate the distance (in pixels) between the calculated and expected
    distance = [sqrt((xc-xe)**2 + (yc-ye)**2) for xc, yc, xe, ye in zip(xcal, ycal, xexp, yexp)]
    d_array = numpy.array(distance)
    
    # Print some stats
    print "Num Reflections: {0}", len(d_array)
    print "Distance (in pixels):"
    print "Min: {0:.2}; Max: {1:.2}; Mean: {2:.2}; Sdev: {3:.2}".format(
        numpy.min(d_array), 
        numpy.max(d_array), 
        numpy.mean(d_array), 
        numpy.std(d_array))

    
def exercise_2(gxparam_path, integrate_path):
    """Execute the second cctbx exercise.
    
    In this exercise, we'll be using the XDS files, with the GXPARM file 
    specifying the geometry of the experiment, and the INTEGRATE file
    specifying the reflection parameters. 
    
    Taking the HKL miller indices and zcal from INTEGRATE.HKL, we will be able
    to calculate the xcal and ycal, the calculated reflection positions. The
    numbers we get out won't be exactly the same because XDS does some
    refinement, but they should be within a pixel or so.
    
    :param gxparam_path: The path to the GXPARM file
    :param integrate_path: The path to the INTEGRATE file
    
    """
    # Read the gxparam file
    gxparm_handle = GxParmFile()
    gxparm_handle.read_file(gxparam_path)

    # Read the integrate file
    integrate_handle = IntegrateFile()
    integrate_handle.read_file(integrate_path)

    # Get the parameters we need
    miller_indices     = integrate_handle.hkl
    xyzcal             = integrate_handle.xyzcal
    wavelength         = gxparm_handle.wavelength
    pixel_size         = gxparm_handle.pixel_size
    detector_origin    = gxparm_handle.detector_origin
    detector_distance  = gxparm_handle.detector_distance
    detector_normal    = gxparm_handle.detector_normal
    beam_direction     = gxparm_handle.beam_vector
    rotation_axis      = gxparm_handle.rotation_axis
    starting_angle     = gxparm_handle.starting_angle
    starting_frame     = gxparm_handle.starting_frame
    oscillation_rangle = gxparm_handle.oscillation_range
    beam_vector        = gxparm_handle.beam_vector
    #unit_cell_a_axis   = gxparm_handle.unit_cell_a_axis
    #unit_cell_b_axis   = gxparm_handle.unit_cell_b_axis
    #unit_cell_c_axis   = gxparm_handle.unit_cell_c_axis
    detector_x_axis    = gxparm_handle.detector_x_axis
    detector_y_axis    = gxparm_handle.detector_y_axis
    
    # Put the detector x and y axes into pixel coordinates
    detector_x_axis = matrix.col(detector_x_axis) / pixel_size[0]
    detector_y_axis = matrix.col(detector_y_axis) / pixel_size[1]
    #detector_origin = detector_distance * matrix.col(detector_normal)
    print wavelength
    # The unit cell parameters have been refined so need to use data in
    # the integrate.hkl file rather than gxparm.xds file.    
    unit_cell_a_axis = (33.914, -36.200, -14.128)
    unit_cell_b_axis = (38.856,  31.328,  13.000)
    unit_cell_c_axis = (-0.543, -19.192,  47.871)

    #iub = matrix.sqr(unit_cell_a_axis + unit_cell_b_axis + unit_cell_c_axis)
    #print iub.inverse()

    xcal = [xyz[0] for xyz in xyzcal]
    ycal = [xyz[1] for xyz in xyzcal]  
    zcal = [xyz[2] for xyz in xyzcal]
   
    # Calculate the current angle of rotation
    rotation_angle = []
    for z in zcal:
        current_frame = z
        rotation_angle.append(starting_angle + oscillation_rangle * (
            current_frame - starting_frame))
   

   
    # Predict the detector coordinates of the given reflections.
    coords = calculate_reflection_detector_coordinates(
        miller_indices, 
        beam_vector, 
        rotation_axis, 
        unit_cell_a_axis, 
        unit_cell_b_axis, 
        unit_cell_c_axis, 
        detector_origin[0], 
        detector_origin[1], 
        detector_distance, 
        detector_x_axis, 
        detector_y_axis, 
        detector_normal, 
        rotation_angle)
        
    # Check the calculated coordinates with those in the file
    xcoords = [c[0] for c in coords]
    ycoords = [c[1] for c in coords]
    validate_calculated_coordinates(xcoords, ycoords, xcal, ycal)
    
    
if __name__ == '__main__':
    
    # Set the XDS files to use
    gxparam_path = '../data/GXPARM.XDS'
    integrate_path = '../data/INTEGRATE.HKL' 
    
    # Execute the second exercise.
    exercise_2(gxparam_path, integrate_path)
    
