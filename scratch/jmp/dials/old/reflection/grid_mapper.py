

class GridMapper:
    """A class to map reflections from the detector to their own specific
    coordinate system grid.

    """

    def __init__(self, gonio, detector, beam, grid, dcs, dcs_to_lcs, image, dp=(4, 4, 1)):
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
        self.dcs = dcs
        self.image = image
        self.dcs_to_lcs = dcs_to_lcs

    def calculate_mask_array(self, rmask, rcs, x, y, z):
        """Calculate the masked array from the reflection mask.

        Args:
            rmask: The reflection mask
            rcs: The reflection coordinate system
            x: The x pixel coordinate
            y: The y pixel coordinate
            z: The z pixel coordinate

        Returns:
            The masked array

        """
        from numpy import zeros, int32
        from scitbx import matrix

        # Calculate the range of detector indices and frames to look at
        x0 = int(x) - self.dp[0]
        x1 = int(x) + self.dp[0]
        y0 = int(y) - self.dp[1]
        y1 = int(y) + self.dp[1]
        z0 = int(z) - self.dp[2]
        z1 = int(z) + self.dp[2]

        # Allocate the mask array
        mask = zeros(shape=(2 * self.dp[2] + 1,
                            2 * self.dp[1] + 1,
                            2 * self.dp[0] + 1), dtype=int32)

        # Loop through all the detector pixels
        for k in range(z0, z1 + 1):
            for j in range(y0, y1 + 1):
                for i in range(x0, x1 + 1):

                    # Calculate the beam vector corresponding to the
                    # detector coordinate.
                    s_dash = matrix.col(self.dcs_to_lcs.apply((i, j)))
#                    s_dash = self.dcs.to_laboratory_coordinate(i, j)
                    s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
                    phi_dash = self.gonio.get_angle_from_frame(k)

                    # Calculate the reflection coordinate from the beam
                    # vector and angle
                    c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash,
                                                                 phi_dash)

                    # Set the mask value
                    mask[k-z0, j-y0, i-x0] = rmask.is_in_mask(c1, c2, c3)

        # return the mask
        return mask

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
        from operator import mul
        import numpy

        # If no maximum is set, then set it to 0.1 N pixels
        if max_iter == None:
            max_iter = int(0.1 * reduce(mul, pixels.shape))

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

    def calculate_e3_fraction(self, rcs, lcs_to_rcs, z):
        """Calculate the fraction of counts contributed by each data frame, j,
        around the reflection to each grid point, v3 in the profile frame.

        First we find and integrate over the set of phi angles covered by each
        data frame around the reflection to get Ij. Then we find the set of phi
        angles covered by each grid coordinate. We then integrate over the
        intersection of the phi angles covered by each data frame and each
        grid point to get Iv3j. The fraction of the counts is then calculated as
        Iv3j / Ij.

        Further details of the method can be found in Kabsch 2010.

        Args:
            rcs: The reflection coordinate system
            z: The z coordinate of the reflection

        Returns:
            An array containing the count fractions. The fraction of counts
            given by frame j to grid coordinate v3 can be found in the array
            by fv3j[v3-v30, j-j0]

        """
        from numpy import zeros, float32
        from math import sqrt, erf

        # Create an array to contain the intensity fractions
        fv3j = zeros(shape=(self.grid.grid_size[2], 2 * self.dp[2] + 1),
            dtype=float32)

        # Get some parameters
        delta3 = self.grid.step_size[2]
        phi0 = self.gonio.starting_angle
        dphi = self.gonio.oscillation_range
        zeta = rcs.zeta
        phi = self.gonio.get_angle_from_frame(z)

        # The range of data frames and grid points to iterate over
        j0 = int(z) - self.dp[2]
        j1 = int(z) + self.dp[2]
        v30 = self.grid.grid_origin[2] - (self.grid.grid_size[2] - 1)
        v31 = (self.grid.grid_size[2] - 1)- self.grid.grid_origin[2]

        # A constant used in the solution to the integrals below.
        sigr2 = 1.0 / (sqrt(2) * (self.grid.sigma_mosaicity / abs(rcs.zeta)))

        # Loop over all j data frames in the region around the reflection
        for j in range(j0, j1 + 1):

            # The data frame j covers the range of phi such that
            # rj = {phi':phi0 + (j-1)dphi <= phi' >= phi0 + jdpi}
            # Therefore the range of phi for j is given as follows.
            aj = phi0 + (j - 1) * dphi
            bj = phi0 + j * dphi

            # Calculate the integral over rj (leaving out scaling factors):
            # I[exp(-(phi' - phi)^2 / (2 sigma^2)]
            Ij = (erf((bj - phi) * sigr2) - erf((aj - phi) * sigr2))

            # Loop over all v3 in the profile grid
            for v3 in range(v30, v31 + 1):

                # The grid coordinate v3 cover the range phi such that
                # rv3 = {phi':(v3 - 0.5)d3 <= (phi' - phi)zeta <= (v3 + 0.5)d3}
                # Therefore the range of phi for v3 is given as follows.
                bv3 = ((v3 - 0.5) * delta3) / zeta + phi
                av3 = ((v3 + 0.5) * delta3) / zeta + phi
                av3, bv3 = (min([av3, bv3]), max([av3, bv3]))

                # We need to integrate over the intersection of sets rv3 and rj
                av3j = max([av3, aj])
                bv3j = min([bv3, bj])

                print Ij, av3j, bv3j

                # If there is no intersection then set the fraction of the
                # counts contributed by fata frame j to grid coordinate v3 to
                # zero, otherwise calculate it as the ratio of the integral
                # over the intersection or rv3 and rj to the integral over rj
                if (av3j > bv3j):
                    fv3j[v3 - v30, j - j0] = 0
                else:
                    Iv3j = erf((bv3j - phi) * sigr2) - erf((av3j - phi) * sigr2)
                    fv3j[v3 - v30, j - j0] = Iv3j / Ij

        # Return the intensity fractions
        return fv3j

    def calculate_e1e2_fraction(self, rcs, lcs_to_rcs, x, y, z, x0, x1, y0, y1, n_div=5):
        """Calculate the fraction of each x, y pixel value that contributes
        to each i, j grid coordinate.

        Loop through all the x, y coordinates and device them into
        n_div x n_div equal sub divisions. Calculate the grid coordinate of
        each sub division and add set the fraction of counts contributed.

        Args:
            rcs: The reflection coordinate system
            x: The x coordinate of the reflection
            y: The y coordinate of the reflection
            z: The z coordinate of the reflection
            x0: The detector x starting coordinate
            x1: The detector x ending coordinate
            y0: The detector y starting coordinate
            y1: The detector y ending coordinate
            n_div: The number of sub-divisions in either direction

        Returns:
            The fraction of counts contributed by each pixel to each grid point.

        """
        from numpy import zeros, arange, float32
        from scitbx import matrix

        # Create an array to contain the intensity fractions
        fv12ij = zeros(shape=(self.grid.grid_size[0],
                              self.grid.grid_size[1],
                              x1 - x0 + 1,
                              y1 - y0 + 1), dtype=float32)

        # Divide counts equally amount subdivisions
        div_fraction = float(n_div * n_div)

        # Loop through all the x, y coordinates
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):

                # Create 5x5 equal sub divisions of the pixel
                for yy in range(0, n_div):
                    for xx in range(0, n_div):

                        # Calculate sub-pixel coordinate
                        xxx = x + xx / float(n_div)
                        yyy = y + yy / float(n_div)

                        # Calculate the beam vector corresponding to the
                        # detector coordinate.
                        s_dash = matrix.col(self.dcs_to_lcs.apply((xxx, yyy)))
                        s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
                        phi_dash = self.gonio.get_angle_from_frame(z)

                        # Calculate the reflection coordinate from the beam
                        # vector and angle
                        c1, c2, c3 = lcs_to_rcs.apply(s_dash, phi_dash)

                        # Find the grid coordinate containing the point
                        gi = self.grid.grid_origin[0] + c1 / self.grid.step_size[0]
                        gj = self.grid.grid_origin[1] + c2 / self.grid.step_size[1]
                        if (0 <= gi < 9 and 0 <= gj < 9):
                            fv12ij[gj,gi,y-y0,x-x0] += 1.0 / div_fraction

        # Return the fractions
        return fv12ij

    def calculate_grid(self, rcs, lcs_to_rcs, x, y, z, x0, x1, y0, y1, z0, z1,
                       reflection_image, fv3j, n_div=5):
        from numpy import zeros, arange, float32
        from scitbx import matrix

        #c1old, c2old, c3old = (0, 0, 0)
        #for x in range(x0, x1 + 1):
        #    s_dash = matrix.col(self.dcs_to_lcs.apply((x, y0)))
        #    s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
        #    phi_dash = self.gonio.get_angle_from_frame(z0)#

            # Calculate the reflection coordinate from the beam
            # vector and angle
        #    c1, c2, c3 = lcs_to_rcs.apply(s_dash, phi_dash)
        #
        #    print (matrix.col((c1, c2, c3)) - matrix.col((c1old, c2old, c3old))).length()
        #    c1old, c2old, c3old = (c1, c2, c3)

        #print 1/0


        # Divide counts equally amount subdivisions
        div_fraction_r = 1.0 / float(n_div * n_div)

        # Loop through all the x, y coordinates
        for y in range(y0, y1 + 1):
            for x in range(x0, x1 + 1):
                for z in range(z0, z1 + 1):

                    # Get the image value
                    value = reflection_image[z-z0, y-y0, x-x0]

                    # Create 5x5 equal sub divisions of the pixel
                    for yy in range(0, n_div):
                        for xx in range(0, n_div):

                            # Calculate sub-pixel coordinate
                            xxx = x + xx / float(n_div)
                            yyy = y + yy / float(n_div)

                            # Calculate the beam vector corresponding to the
                            # detector coordinate.
                            s_dash = matrix.col(self.dcs_to_lcs.apply((xxx, yyy)))
                            s_dash = (1.0 / self.beam.wavelength) * s_dash.normalize()
                            phi_dash = self.gonio.get_angle_from_frame(z)

                            # Calculate the reflection coordinate from the beam
                            # vector and angle
                            c1, c2, c3 = lcs_to_rcs.apply(s_dash, phi_dash)

                            # Find the grid coordinate containing the point

                            gi = self.grid.grid_origin[0] + c1 / self.grid.step_size[0]
                            gj = self.grid.grid_origin[1] + c2 / self.grid.step_size[1]

                            for gk in range(0, 9):

                                if (0 <= gi < 9 and 0 <= gj < 9 and 0 <= gk < 9):
                                    fraction = fv3j[z-z0,gk] * div_fraction_r
                                    self.grid.grid_data[0,gk,gj,gi] += value * fraction


    def map_reflection(self, index, rcs, lcs_to_rcs, rmask, x, y, z, s1, phi, wavelength):
        """Map the reflection from the detector to the profile grid.

        Loop through all the pixels around the reflection in the detector
        image volume. Transform the detector coordinate to the reflection
        profile coordinate system, then distribute the value of the pixel
        to the reflection grid.

        Args:
            index: The reflection number
            rcs: The reflection coordinate system
            rmask: The reflection mask
            x: The reflection x detector coordinate
            y: The reflection y detector coordinate
            z: The reflection z detector coordinate (frame number)

        """
        from scitbx.array_family import flex
        from scitbx import matrix
        from numpy import where, arange, logical_not, array
        from math import sqrt

        if (not (self.dp[0] <= x < self.image.shape[2] - self.dp[0] - 1) or
            not (self.dp[1] <= y < self.image.shape[1] - self.dp[1] - 1) or
            not (self.dp[2] <= z - self.gonio.starting_frame < self.image.shape[0] - self.dp[2] - 1)):
            return

        # Calculate the range of detector indices and frames to look at
        x0 = int(x) - self.dp[0]
        x1 = int(x) + self.dp[0]
        y0 = int(y) - self.dp[1]
        y1 = int(y) + self.dp[1]
        z0 = int(z) - int(self.gonio.starting_frame) - self.dp[2]
        z1 = int(z) - int(self.gonio.starting_frame) + self.dp[2]

        # Get the reflection image
        reflection_image = self.image[z0:z1+1,y0:y1+1,x0:x1+1]

        # Calculate the reflection mask
        #reflection_mask = self.calculate_mask_array(rmask, rcs, x, y, z)

        # Calculate the background intensity and subtract from image. Ensure
        # that counts are not less than zero
        background = self.calculated_background_intensity(reflection_image)

        #background = self.calculated_background_intensity(
        #    reflection_image[where(reflection_mask != 0)])
        reflection_image -= background
        reflection_image[where(reflection_image < 0.0)] = 0.0
        #reflection_image[where(reflection_mask == 0)] = 0.0


        # Calculate the fraction of counts constributed by each data frame
        # around the reflection to each v3 grid point in the reflection grid
        #fv3j = self.calculate_e3_fraction(rcs, lcs_to_rcs, z)
        from dials_jmp.geometry.algorithm import XdsTransform
        from dials_jmp.geometry.transform import FromDetectorToXds
        from dials_jmp.geometry.algorithm import XdsTransformGrid

        dcs_to_xcs = FromDetectorToXds(self.dcs.in_si_units(self.detector.pixel_size),
                                       self.detector.origin,
                                       self.detector.distance,
                                       rcs,
                                       s1,
                                       phi,
                                       wavelength)

        grid_n = (int(self.grid.grid_size[0]/2), int(self.grid.grid_size[1]/2),int(self.grid.grid_size[2]/2))

        #frame_angles = flex.double(array(
        #    map(lambda z: self.gonio.get_angle_from_frame(z),
        #        arange(0, 180))))

        xds_grid = XdsTransformGrid(1, (4, 4, 4),
                                    self.grid.sigma_divergence,
                                    self.grid.sigma_mosaicity)

        print xds_grid.n_reflections
        print xds_grid.size
        print xds_grid.origin
        print xds_grid.step_size
        print xds_grid.sigma_mosaicity
        print xds_grid.delta_mosaicity


        xds_trans = XdsTransform(xds_grid,
                                 flex.int(self.image),
                                 self.image.shape[::-1],
                                 self.detector,
                                 self.beam,
                                 self.gonio,
                                 self.dp)

        import time

        st = time.time()
#        grid = xds_trans.calculate3(lcs_to_rcs, (x, y, z), phi, rcs.zeta)
        xds_trans.calculate(0, (x, y, z), s1, phi)

        grid = xds_grid.data

        print "Time Taken: ", time.time() - st


        grid = grid.as_numpy_array()

        grid.shape = (1, 9, 9, 9)
        self.grid.grid_data = grid
#        self.grid.grid_data[index,:,:,:] = grid
        #print grid
        #print 1/0
        #fv3j = xds_trans.calculate(z, phi, rcs.zeta)
        #fv3j = fv3j.as_numpy_array()
        #fv3j.shape = (2*self.dp[2]+1, self.grid.grid_size[2])

        #fv1v2xy = self.calculate_grid(rcs, lcs_to_rcs, x, y, z, x0, x1, y0, y1, z0, z1,
        #                                reflection_image, fv3j)
        #print 1/0
        # Calculate the fraction of counts from each pixel coordinate that
        # goes to each grid coordinate
        #fv1v2xy = self.calculate_e1e2_fraction(rcs, lcs_to_rcs, x, y, z, x0, x1, y0, y1, n_div=5)

        # Loop through all the detector coordinates.
        # TODO this needs to be fixed (6 nested loops!)
        #for z in range(z0, z1 + 1):
        #    for y in range(y0, y1 + 1):
        #        for x in range(x0, x1 + 1):
        #
        #            # Get the image value
        #            value = reflection_image[z-z0, y-y0, x-x0]##

                    # Loop through all the grid coordinates
        #            for k in range(0, self.grid.grid_size[2]):
        #                for j in range(0, self.grid.grid_size[1]):
        #                    for i in range(0, self.grid.grid_size[0]):

                                # Calculate the pixel fraction and add the
        #                        # image counts
        #                        fraction = fv3j[z-z0,k] * fv1v2xy[j,i,y-y0,x-x0]
        #                        self.grid.grid_data[index,k,j,i] += value * fraction
