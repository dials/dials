
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
    
    def is_data_normally_distributed(self, data, delta=0.1):
        """Check if the data is normally distributed.
        
        Calculate the percentage number of data values within +- 1, 2 and 3 
        standard deviations of the mean value. Values that are normally 
        distributed will have:
            ~68.2% of points within +- 1 standard deviation
            ~95.4% of points within +- 2 standard deviation
            ~99.7% of points within +- 3 standard deviation                    

        If the percentages match (within +- delta) for the given data then the
        function returns true, otherwise returns false.

        Args:
            data: The array of data
            delta: The allowed difference
            
        Returns:
            True/False

        """    
        import numpy
    
        # Calculate the mean and standard deviation of the data
        n_data = len(data)
        mean = numpy.mean(data)
        sdev = numpy.std(data)
        
        # Calculate the mean +- [1|2|3] standard deviations
        sd1p, sd1m = (mean + 1 * sdev, mean - 1 * sdev)
        sd2p, sd2m = (mean + 2 * sdev, mean - 2 * sdev)
        sd3p, sd3m = (mean + 3 * sdev, mean - 3 * sdev)

        # Calculate percentage of data within [1|2|3] standard deviations
        perc1 = len(numpy.where((sd1m <= data) & (data < sd1p))[0]) / n_data
        perc2 = len(numpy.where((sd2m <= data) & (data < sd2p))[0]) / n_data
        perc3 = len(numpy.where((sd3m <= data) & (data < sd3p))[0]) / n_data
                
        # Return whether the data is normally distributed
        return ((0.682 - delta <= perc1 <= 0.682 + delta) and
                (0.954 - delta <= perc2 <= 0.954 + delta) and
                (0.997 - delta <= perc2 <= 0.997 + delta))
    
    def calculated_background_intensity(self, pixels, delta=0.1, 
                                        max_iter=None):
        """Calculate the background intensity.
        
        Sort the pixels in order of ascending intensity. Then check if the
        intensities are normally distributed. If not then remove the pixel
        with the highest intensity from the list and check again. Keep going
        untill the list of pixels is normally distributed, or the maximum 
        number of iterations is reached. Return the mean of the values as the
        background intensity.
        
        Args:
            pixels: The list of pixels
            delta: The allowed difference from normal
            max_iter: The maximum number of iterations
           
        Returns:
            The background intensity value
        
        """        
        import numpy

        # If no maximum is set, then set it to 0.1 N pixels
        if max_iter == None:
            max_iter = int(0.1 * len(pixels.flat))

        # Sort the pixels into ascending intensity order    
        sorted_pixels = pixels.flat[:]
        sorted_pixels.sort()

        # Check if the data is normally distributed. If it is not, then remove
        # a value of high intensity and keep looping until it is. If the number
        # of iterations exceeds the maximum then exit the loop.
        num_iter = 1
        while not self.is_data_normally_distributed(sorted_pixels, delta=delta):
            sorted_pixels = sorted_pixels[:-1]
            if num_iter >= max_iter:
                break
            num_iter += 1
                     
        # Return the mean of the remaining pixels as the background intensity
        return numpy.mean(sorted_pixels)      
    
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
        
        # Calculate the background intensity
        background = self.calculated_background_intensity(
            self.detector.image[int(z)-1:int(z)+1+1, 
                                int(y)-4:int(y)+4+1, 
                                int(x)-4:int(x)+4+1])
        
        # Get the detector coordinate system
        dcs = self.detector.coordinate_system

        # Calculate the range of detector indices and frames to look at
        x0 = int(x) - self.dp[0]
        x1 = int(x) + self.dp[0]
        y0 = int(y) - self.dp[1]
        y1 = int(y) + self.dp[1]
        z0 = int(z) - self.dp[2]
        z1 = int(z) + self.dp[2]
        
        import numpy
        from scipy.special import erf
        from math import sqrt
        #all_phi_dash = map(self.gonio.get_angle_from_frame, 
        #                   self.gonio.z0 + 
        #                   numpy.arange(0, 
        #                        20*self.detector.image.shape[0])/20)
        
        
        
        
        delta3 = self.grid.step_size[2]
        #print delta3
        zeta = rcs.zeta
        phi = rcs.phi
        phi0 = self.gonio.phi0
        dphi = self.gonio.dphi
        
        sigmar2 = sqrt(2) * (self.grid.sigma_mosaicity / abs(rcs.zeta)) 
        
        fv3j = numpy.zeros(shape=(9,3), dtype=numpy.float32)
        for v3 in range(-4, 4+1):
            for j in range(-1, 1+1):

                av3 = ((v3 - 0.5) * delta3) / zeta + phi
                bv3 = ((v3 + 0.5) * delta3) / zeta + phi
                av3, bv3 = (min([av3, bv3]), max([av3, bv3]))
                
                aj = phi0 + (int(z) + j - 1) * dphi
                bj = phi0 + (int(z) + j) * dphi
                
                av3j = max([av3, aj])
                bv3j = min([bv3, bj])

                if (v3 == -1):
                    print av3, bv3, aj, bj, av3j, bv3j

                num = (erf((bv3j - phi) / sigmar2) - erf((av3j - phi) / sigmar2))
                den = (erf((bj - phi) / sigmar2) - erf((aj - phi) / sigmar2))

                #print aj, bj, phi, erf((bj - phi) / sigmar2), erf((aj - phi) / sigmar2)
                #print num, den                
                if (av3j > bv3j):
                    fv3j[v3+4,j+1] = 0
                else:
                    fv3j[v3+4,j+1] = num / den
        
        print fv3j
        
        from scitbx import matrix
        
        for z in range(4, 6 + 1):
            for y in range(y0, y1 + 1):
                for x in range(x0, x1 + 1):
                    
                    for yy in numpy.arange(0, 10) / 10.0:
                        for xx in numpy.arange(0, 10) / 10.0:
                            
                            xxx = x + xx
                            yyy = y + yy
                    
                            # Calculate the beam vector corresponding to the detector 
                            # coordinate.
                            s_dash = matrix.col(dcs.to_laboratory_coordinate(xxx, yyy))
                            s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
                            phi_dash = self.gonio.get_angle_from_frame(z)
                            
                            
                            # Calculate the reflection coordinate from the beam vector
                            c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash, 
                                                                         phi_dash)
                                                                         

                            for gk in range(0,9):
                                gi = self.grid.grid_origin[0] + c1 / self.grid.step_size[0]
                                gj = self.grid.grid_origin[1] + c2 / self.grid.step_size[1]
                                #gk = self.grid.grid_origin[2] + c3 / self.grid.step_size[2]
                                c3 = (gk - self.grid.grid_origin[2]) * self.grid.step_size[2]
                                if (0 <= gi < 9 and 0 <= gj < 9 and 0 <= gk < 9):
                                    value = (self.detector.image[z-1, y, x] - background)
                                    if (value < 0):
                                        value = 0
                                    value = (self.detector.image[z-1, y, x] - background) * fv3j[gk,z-4] / 100.0
                                    self.grid.grid_data[0,gk,gj,gi] += value
        
        
                                    # Add the count from the pixel to the reflection grid
                                    #value = self.detector.image[z-z0, y, x]
                                    #self.grid.add_point_count(index, (c1, c2, c3), value)
                    
#        print numpy.array(fv3j * 100000, dtype=numpy.int32)
        
        
        #set_v3 = {}
        #for phi_dash in all_phi_dash:
        #    set_v3[phi_dash] = []
        #    for v3 in range(-int(self.grid.grid_size[2]/2), int(self.grid.grid_size[2]/2)+1):
        #         if (v3 - 0.5) * delta3 <= (phi_dash - phi) * zeta <= (v3 + 0.5) * delta3:
        #            set_v3[phi_dash].append(v3)
        
        #ji = range(1, 181)
        #set_j = {}
        #for phi_dash in all_phi_dash:
        #    set_j[phi_dash] = []
        #    for j in ji:
        #        if (1 + (j - 1) * 1.0 <= phi_dash <= 1 + j * 1.0):
        #            set_j[phi_dash].append(j)
        
        #new_set = []
        #for phi_dash in all_phi_dash:
        #    if (set_j[phi_dash] and set_v3[phi_dash]):
        #        new_set.append(phi_dash)
        
        #print set_j
        #print set_v3
        #new_set_j = []
        #for phi_dash in all_phi_dash:
        #    if (set_j[phi_dash]):
        #        new_set_j.append(phi_dash)
        
        #new_set_v3 = []
        #for phi_dash in all_phi_dash:
        #    if (set_v3[phi_dash]):
        #        new_set_v3.append(phi_dash)
        
        
        #print "Set_v3", new_set_v3
        #print "Intersection of set_j and set_v3", new_set
        #for p in new_set:
        #    print "Phi_dash {0}, set_j[phi_dash] {1}, set_v3[phi_dash] {2}".format(
        #        p, set_j[p], set_v3[p])
        
        #print "e3[phi_dash = 5] {0}, e3[phi_dash = 6]{1}".format(
        #    rcs.from_diffracted_beam_vector(rcs.s1, 5)[2], 
        #    rcs.from_diffracted_beam_vector(rcs.s1, 6)[2])
        
        #print "e3[grid] = {0}".format(
        #    map(lambda x: self.grid.calculate_grid_coordinate((0, 0, x))[2], 
        #        range(0,9)))
        
        
        
        #for phi_dash in all_phi_dash:
        #    c1, c2, c3 = rcs.from_diffracted_beam_vector(rcs.s1, phi_dash)
        #    print c3
        
        #print 5.5*self.grid.sigma_mosaicity * (180 / (3.14159))
       # 
        #for v3 in range(-int(self.grid.grid_size[2]/2), int(self.grid.grid_size[2]/2)+1):
        #    #v1, v2, v3 = self.grid.calculate_grid_coordinate((0, 0, k))
        #    self.calculate_fraction_of_counts(0, v3, rcs, self.grid, all_phi_dash)
    
        
        # Loop through all the detector pixels in the selected imaqes
        #for z in range(z0, z1 + 1):
        #    for y in range(y0, y1 + 1):
        #        for x in range(x0, x1 + 1):
        #            
                    # Calculate the beam vector corresponding to the detector 
                    # coordinate.
        #            s_dash = dcs.to_laboratory_coordinate(x, y)
        #            s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
        #            phi_dash = self.gonio.get_angle_from_frame(z)
                    
                    # Calculate the reflection coordinate from the beam vector
        #            c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash, 
        #                                                         phi_dash)
        
                    # Add the count from the pixel to the reflection grid
        #            value = self.detector.image[z-z0, y, x]
        #            self.grid.add_point_count(index, (c1, c2, c3), value)
                    

