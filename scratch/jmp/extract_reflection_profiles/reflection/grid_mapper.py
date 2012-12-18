
class GridMapper:
    """A class to map reflections from the detector to their own specific 
    coordinate system grid.
        
    """
    
    def __init__(self, gonio, detector, beam, grid, dp=(4, 4, 4)):
        """Set all the parameters needed.
        
        Args:
            gonio: The goniometer
            image: The 3D image volume
            dcs: The detector coordinate system
            grid: The reflection profile grid
            wavelength: The wavelength of the x-rays
            dp: (dx, dy, dz) The number of pixels either side to use.
        
        """
        self.gonio = gonio
        self.detector = detector
        self.beam = beam
        self.grid = grid
        self.dp = dp
    
    def calculate_fraction_of_counts(self, j, v3, rcs, grid, all_phi_dash):

        delta3 = grid.step_size[2]
        zeta = rcs.zeta
        phi = rcs.phi
        
        sigma2 = 2.0 * (grid.sigma_mosaicity / abs(rcs.zeta))**2
        
        for phi_dash in all_phi_dash:
            print (v3 - 0.5) * delta3 <= (phi_dash - phi) * zeta <= (v3 + 0.5) * delta3
        
    
    def map_reflection(self, index, rcs, x, y, z):
        """Map the reflection from the detector to the profile grid.
        
        Loop through all the pixels around the reflection in the detector
        image volume. Transform the detector coordinate to the reflection
        profile coordinate system, then distribute the value of the pixel
        to the reflection grid.
        
        Args:
            index: The reflection number
            rcs: The reflection coordinate system
            x: The reflection x detector coordinate
            y: The reflection y detector coordinate
            z: The reflection z detector coordinate (frame number)
        
        """
        # Get the detector coordinate system
        dcs = self.detector.coordinate_system

        # Calculate the range of detector indices and frames to look at
        x0 = int(x) - self.dp[0]
        x1 = int(x) + self.dp[0]
        y0 = int(y) - self.dp[1]
        y1 = int(y) + self.dp[1]
        z0 = int(z) - self.dp[2]
        z1 = int(z) + self.dp[2]
        
        # Loop through all the detector pixels in the selected imaqes
        for z in range(z0, z1 + 1):
            for y in range(y0, y1 + 1):
                for x in range(x0, x1 + 1):
                    
                    # Calculate the beam vector corresponding to the detector 
                    # coordinate.
                    s_dash = dcs.to_laboratory_coordinate(x, y)
                    s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
                    phi_dash = self.gonio.get_angle_from_frame(z)
                    
                    # Calculate the reflection coordinate from the beam vector
                    c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash, 
                                                                 phi_dash)
        
                    
                    import numpy
                    all_phi_dash = map(self.gonio.get_angle_from_frame, 
                                       self.gonio.z0 + 
                                       numpy.arange(0, 
                                            self.detector.image.shape[0]))
                    
                    self.calculate_fraction_of_counts(0, 0, rcs, self.grid, all_phi_dash)
    
                    # Add the count from the pixel to the reflection grid
                    value = self.detector.image[z-z0, y, x]
                    self.grid.add_point_count(index, (c1, c2, c3), value)
                    

