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
        pass
    
    def read_file(self, filename):
        """Read the INTEGRATE.HKL file.
        
        See http://xds.mpimf-heidelberg.mpg.de/html_doc/xds_files.html for more
        information about the file format.
        
        :param filename: The path to the file
        
        """
        # Read the lines from the file
        lines = open(filename, 'r').readlines()
        
        # Parse each line of the file. Store each column in the file in a list
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
        for l in lines:
            
            # If line begins with a '!' then ignore
            if l.lstrip().startswith('!'):
                continue
        
            # Split the tokens
            tokens = l.split()
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
                                              y0, f, d1, d2, d3, phi0, dphi, ub):

   
    m2 = matrix.col(m2)
    s0 = matrix.col(s0)
    b1 = matrix.col(b1)
    b2 = matrix.col(b2)
    b3 = matrix.col(b3)
    d1 = matrix.col(d1)
    d2 = matrix.col(d2)
    d3 = matrix.col(d3)

    b1_star = b2.cross(b3) / (b1.dot((b2.cross(b3))))
    b2_star = b3.cross(b1) / (b1.dot((b2.cross(b3))))
    b3_star = b1.cross(b2) / (b1.dot((b2.cross(b3))))
    
       
    #print indices[30513]
   
    coords = []
    for i in [1]:
        h, k, l = indices[i]
        phi = phi0[i]

        #r = m2.axis_and_angle_as_r3_rotation_matrix(phi, deg=True)

        p_star0 = h*b1_star + k*b2_star + l*b3_star
        p_star = p_star0.rotate(axis=m2, angle=phi, deg=True)
        #p_star = r * ub * matrix.col((h, k, l))
        s = s0 + p_star

        print "H, K, L:", h, k, l
        print "Phi:", phi
        #print "ub:", ub
        print "m2:", tuple(m2)
        print "b1", tuple(b1)
        print "b2", tuple(b2)
        print "b3", tuple(b3)
        print "b1_star", tuple(b1_star)
        print "b2_star", tuple(b2_star)
        print "b3_star", tuple(b3_star)
        print "p_star0:", tuple(p_star0)
        print "p_star:", tuple(p_star)
        print "S0:", tuple(s0)
        print "S:", tuple(s)
        print "Length S0:", s0.length()
        print "Length S:", s.length()
       
       
#        s_dot_d3 = s.dot(d3)
#        print "S.d3:", s_dot_d3
#        if (s_dot_d3 > 0):
#            x = x0 + f*s.dot(d1) / s_dot_d3
#            y = y0 + f*s.dot(d2) / s_dot_d3
#            print "X, Y:", x, y
#            coords.append((x, y))

#    print coords[0]

def validate_calculated_coordinates():
    pass

    
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
    detector_origin    = gxparm_handle.detector_origin
    detector_distance  = gxparm_handle.detector_distance
    detector_normal    = gxparm_handle.detector_normal
    beam_direction     = gxparm_handle.beam_vector
    rotation_axis      = gxparm_handle.rotation_axis
    starting_angle     = gxparm_handle.starting_angle
    starting_frame     = gxparm_handle.starting_frame
    oscillation_rangle = gxparm_handle.oscillation_range
    beam_vector        = gxparm_handle.beam_vector
    crystal_x_axis     = gxparm_handle.unit_cell_a_axis
    crystal_y_axis     = gxparm_handle.unit_cell_b_axis
    crystal_z_axis     = gxparm_handle.unit_cell_c_axis
    detector_x_axis    = gxparm_handle.detector_x_axis
    detector_y_axis    = gxparm_handle.detector_y_axis
    
    # Calculate the UB matrix
    iub_matrix = matrix.sqr(crystal_x_axis + crystal_y_axis + crystal_z_axis)
    ub_matrix = iub_matrix.inverse()
    
    zcal = [xyz[2] for xyz in xyzcal]
    
    print xyzcal[30513]
    print "1/Lambda: ", 1.0 / wavelength
   
    # Calculate the current angle of rotation
    rotation_angle = []
    for z in zcal:
        current_frame = z
        rotation_angle.append(starting_angle + oscillation_rangle * (
            current_frame - starting_frame))
        
    
    # Calculate the rotation matrix
    #rotation_matrix = axis_angle_to_matrix(rotation_axis, rotation_angle)
    
    # Predict the detector coordinates of the given reflections.
    coords = calculate_reflection_detector_coordinates(
        miller_indices, 
        beam_vector, 
        rotation_axis, 
        crystal_x_axis, 
        crystal_y_axis, 
        crystal_z_axis, 
        detector_origin[0], 
        detector_origin[1], 
        detector_distance, 
        detector_x_axis, 
        detector_y_axis, 
        detector_normal, 
        rotation_angle, 
        oscillation_rangle,
        ub_matrix)
        
    # Check the calculated coordinates with those in the file
    validate_calculated_coordinates()
    
    
if __name__ == '__main__':
    
    # Set the XDS files to use
    gxparam_path = '../data/GXPARM.XDS'
    integrate_path = '../data/INTEGRATE.HKL' 
    
    # Execute the second exercise.
    exercise_2(gxparam_path, integrate_path)
    