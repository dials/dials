
class Reflection:
    """A class to represent a reflection"""

    def __init__(self, hkl, phi):
        """Initialise the reflection.
        
        Args:
            hkl: The miller indices
            phi: The rotation angle
        
        """
        self.hkl = hkl
        self.phi = phi
        self.s1  = None
        self.xyz = None

    def in_phi_range(self, phi0, phi1):
        """Check whether phi is within the range of phi given.
        
        Args:
            phi0: The starting angle
            phi1: The ending angle
            
        Returns:
            True/False
        
        """
        return (self.phi >= phi0 and self.phi <= phi1)

    def in_detector_volume(self, volume):
        """Check that the detector coordinates of the reflection are within the
        image volume.
        
        Args:
            volume: The image volume ((x0, x1), (y0, y1), (z0, z1))
            
        Returns:
            True/False
        
        """
        return (volume[0][0] <= self.xyz[0] < volume[0][1] and
                volume[1][0] <= self.xyz[1] < volume[1][1] and
                volume[2][0] <= self.xyz[2] < volume[2][1])

    def calculate_detector_coordinates(self, m2, s0, x0, y0, f, d1, d2, d3, 
                                       phi0, dphi, z0, ub):
        """Calculate the pixel coordinates of the reflection.
        
        Calculate the diffracted beam wave vector and find where it intersects
        with the dectector plane.
        
        Args:
            m2: The rotation axis
            s0: The incident beam wave vector
            x0: The x coordinate of the detector origin
            y0: The y coordinate of the detector origin
            f: The distance from the crystal to the detector
            d1: The detector x axis vector
            d2: The detector y axis vector
            d3: The detector normal unit vector
            phi0: The starting angle
            dphi: The oscillation range
            z0: The starting frame
            ub: The UB matrix
        
        Returns:
            True/False
                
        """
        from scitbx import matrix
        
        # Calculate the reciprocal lattice vector and rotate it by phi about
        # the rotation axis. Then calculate the diffracted beam vector, s1.
        rub = m2.axis_and_angle_as_r3_rotation_matrix(self.phi, deg=True) * ub
        p_star = rub * matrix.col(self.hkl)
        self.s1 = s0 + p_star
    
        # Calculate the point (in pixels) at which the diffracted beam vector 
        # intersects the detector plane. (only when s1.d3 > 0)
        s1_dot_d3 = self.s1.dot(d3)
        if (s1_dot_d3 > 0):
            x = x0 + f * self.s1.dot(d1) / s1_dot_d3
            y = y0 + f * self.s1.dot(d2) / s1_dot_d3
            z = (self.phi - phi0) / dphi + z0
            self.xyz = (x, y, z)
            return True

        # Return false
        return False
        
    def calculate_detector_coordinates2(self, gonio, detector, ub_matrix):
        """Calculate the pixel coordinates of the reflection
        
        Args:
            gonio: The goniometer
            detector: The detector
            ub_matrix: The UB matrix
            
        Returns:
            True/False
        
        """
        return self.calculate_detector_coordinates(
            gonio.coordinate_system.m2, gonio.coordinate_system.s0, 
            detector.pixel_origin[0], detector.pixel_origin[1], 
            detector.distance, 
            detector.coordinate_system.d1_1, 
            detector.coordinate_system.d2_1, 
            detector.coordinate_system.d3, 
            gonio.phi0, gonio.dphi, gonio.z0, 
            ub_matrix)
                                       
                                       
                                       
                                       
                                       
