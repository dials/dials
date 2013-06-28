from __future__ import division

class Grid:
    """A class to hold the reflection profiles transformed into the specific
    reflection coordinate system.
    """

    def __init__(self, num_grids, grid_size=None, grid_origin=None,
                 step_size=None, sigma_divergence=0.0, sigma_mosaicity=0.0,
                 n_sigma_divergence=10, n_sigma_mosaicity=10):
        """Initialise the reflection profile grid.

        This class is used to hold grid data in the transformed reflection
        coordinate system. It is convenient to hold all reflection profiles
        in one place rather than allocate small grids for each reflection.

        The step size of the grid can either be set explicitly or can be
        calculated from the beam divergence and the mosaicity.

        Args:
            num_grids: The number of grids (i.e. 1 for each reflection)
            grid_size: The size of the individual grids (defaults to (9, 9, 9))
            grid_origin: The origin of each (defaults to centre)
            step_size: The step size of the grid
            sigma_divergence: The standard deviation of the beam divergence
            sigma_mosaicity: The standard deviation of the mosaicity
            n_sigma_divergence: Number of standard deviation of beam divergence
            n_sigma_mosaicity: Number of standard deviation of mosaicity

        """
        import numpy

        # Set the default values for the grid size and grid origin
        if grid_size == None:
            grid_size = (9, 9, 9)

        if grid_origin == None:
            grid_origin = tuple(map(lambda x: int(x / 2), grid_size))

        # If the step size is not explicitly set, calculate it from the
        # standard deviation of the beam divergence and the mosaicity.
        if step_size == None:

            # If the sigmas are not set, then raise an exception
            if sigma_divergence <= 0.0 or sigma_mosaicity <= 0.0:
                raise "Can't calculate step size"

            # Otherwise, calculate the step size of the grid.
            step_size = self._calculate_step_size(grid_size, sigma_divergence,
                sigma_mosaicity, n_sigma_divergence, n_sigma_mosaicity)

        # Save the step size grid size and number of grids
        self.step_size = step_size
        self.grid_size = grid_size
        self.num_grids = num_grids
        self.grid_origin = grid_origin
        self.sigma_divergence = sigma_divergence
        self.sigma_mosaicity = sigma_mosaicity

        # Allocate the grid data
        grid_shape = (num_grids,) + grid_size[::-1]
        self.grid_data = numpy.zeros(shape=grid_shape, dtype=numpy.float32)

    def estimated_reflection_intensity(self, sigma_divergence, sigma_mosaicity):
        """Calculate the estimated fraction of the total reflection intensity
        held by each element in the profile grids.

        Args:
            sigma_divergence: The standard deviation of the beam divergence
            sigma_mosaicity: The standard deviation of the mosaicity

        Returns:
            A grid of the same dimensions with the fraction of intensity

        """
        from math import exp, sqrt, pi
        import numpy

        # Allocate a grid the same size as the reflection profile grids
        grid_data = numpy.zeros(shape=self.grid_size, dtype=numpy.float32)

        # Set the volume element to the step size
        de1 = self.step_size[0]
        de2 = self.step_size[1]
        de3 = self.step_size[2]

        # Save the denominator of the gaussian function
        rtau_sigma_d = sqrt(2 * pi) * sigma_divergence
        rtau_sigma_m = sqrt(2 * pi) * sigma_mosaicity

        # Save 2*sigma**2
        sigma_d_2 = 2 * sigma_divergence * sigma_divergence
        sigma_m_2 = 2 * sigma_mosaicity * sigma_mosaicity

        # For each element in the grid, calculate the estimated fraction of
        # the total reflection intensity.
        for k in range(0, self.grid_size[2]):
            for j in range(0, self.grid_size[1]):
                for i in range(0, self.grid_size[0]):
                    c1, c2, c3 = self.calculate_grid_coordinate((i, j, k))
                    w1 = de1 * exp(-c1 * c1 / sigma_d_2) / rtau_sigma_d
                    w2 = de2 * exp(-c2 * c2 / sigma_d_2) / rtau_sigma_d
                    w3 = de3 * exp(-c3 * c3 / sigma_m_2) / rtau_sigma_m
                    grid_data[k, j, i] = w1 * w2 * w3

        # Return the grid
        return grid_data

    def _calculate_step_size(self, grid_size, sigma_divergence, sigma_mosaicity,
                             n_sigma_divergence, n_sigma_mosaicity):
        """Calculate the step size from the beam divergence and mosaicity

        Args:
            grid_size: The size of the grid
            sigma_divergence: The standard deviation of the beam divergence
            sigma_mosaicity: The standard deviation of the mosaicity
            n_sigma_divergence: Number of standard deviation of beam divergence
            n_sigma_mosaicity: Number of standard deviation of mosaicity

        Returns:
            The step sizes (de1, de2, de3)

        """
        delta_divergence = n_sigma_divergence * sigma_divergence
        delta_mosaicity = n_sigma_mosaicity * sigma_mosaicity
        return (delta_divergence / grid_size[0],
                delta_divergence / grid_size[1],
                delta_mosaicity  / grid_size[2])

    def get_grid_coordinates(self):
        """Return the e1, e2 and e3 coordinates of the grids

        Returns:
            A tuple of grids containing the grid coordinates

        """
        import numpy

        # Allocate arrays for the grid coordinates
        grid_e1 = numpy.zeros(shape=self.grid_size, dtype=numpy.float32)
        grid_e2 = numpy.zeros(shape=self.grid_size, dtype=numpy.float32)
        grid_e3 = numpy.zeros(shape=self.grid_size, dtype=numpy.float32)

        # Loop through the grid and calculate the coordinates
        for k in range(0, self.grid_size[2]):
            for j in range(0, self.grid_size[1]):
                for i in range(0, self.grid_size[0]):
                    c1, c2, c3 = self.calculate_grid_coordinate((i, j, k))
                    grid_e1[k, j, i] = c1
                    grid_e2[k, j, i] = c2
                    grid_e3[k, j, i] = c3

        # Return the coordinates
        return (grid_e1, grid_e2, grid_e3)

    def get_grid(self, grid_number):
        """Get a reflection grid.

        Args:
            grid_number: The number of grid to return

        Returns:
            The requested grid data.

        """
        return self.grid_data[grid_number]

    def reset(self, grid_number=None):
        """Reset one or all grids to zero

        Args:
            grid_number: The grid to reset

        """
        if grid_number == None:
            self.grid_data[:,:,:,:] = 0
        else:
            self.grid_data[grid_number,:,:,:] = 0

    def normalize(self, grid_number=None):
        """Normalize the grid values to between 0.0 and 1.0

        Args:
            grid_number: The grid to normalize

        """
        import numpy

        if grid_number == None:
            for grid_number in range(0, self.num_grids):
                max_grid = numpy.max(self.grid_data[grid_number,:,:,:])
                self.grid_data[grid_number,:,:,:] /= max_grid
        else:
            max_grid = numpy.max(self.grid_data[grid_number,:,:,:])
            self.grid_data[grid_number,:,:,:] /= max_grid

    def calculate_grid_coordinate(self, index):
        """Calculate the coordinate of a grid point

        Args:
            index: The grid index (i, j, k)

        Returns:
            The grid coordinate (x, y, z)

        """
        return (self.step_size[0] * (index[0] - self.grid_origin[0]),
                self.step_size[1] * (index[1] - self.grid_origin[1]),
                self.step_size[2] * (index[2] - self.grid_origin[2]))

    def _is_point_cube_in_grid(self, index):
        """Check if the cube of grid points specified by i, j, k, i+1, j+1, k+1
        is within the grid dimensions

        Args:
            index: The grid index (i, j, k)

        Returns:
            True/False

        """
        return (0 <= index[0] < self.grid_size[0] - 1 and
                0 <= index[1] < self.grid_size[1] - 1 and
                0 <= index[2] < self.grid_size[2] - 1)

    def _calculate_point_weight(self, point, grid_point):
        """Calculate the weight the give to the grid point. Called by the
        function 'add_point_count'

        Args:
            point: The point at which the count is known
            grid_point: The grid point of a point in the cube around point

        Returns:
            The weight (between 0.0 and 1.0) to give to the grid point

        """
        return (abs(1.0 - abs(grid_point[0] - point[0]) / self.step_size[0]) *
                abs(1.0 - abs(grid_point[1] - point[1]) / self.step_size[1]) *
                abs(1.0 - abs(grid_point[2] - point[2]) / self.step_size[2]))

    def add_point_count(self, grid_number, point, count):
        """Add a point count to the grid

        When a pixel is mapped from the detector to the reflection coordinate
        system, the point is calculated but will probably not coincide
        directly with a grid point. Therefore, the count needs to be distributed
        to the eight surrounding grid points. The least biased way of doing this
        is by maximising the Shannon entropy of the distributed counts.

        Args:
            grid_number: The number of the grid to add the count to
            point: The point at which the count is known
            count: The count

        """
        from math import floor

        # Get the starting indices
        i0 = self.grid_origin[0] + point[0] / self.step_size[0]
        j0 = self.grid_origin[1] + point[1] / self.step_size[1]
        k0 = self.grid_origin[2] + point[2] / self.step_size[2]

        # Ensure cube around point is within grid
        if self._is_point_cube_in_grid((i0, j0, k0)):

            # Now calculate the integer cube indices
            i0, j0, k0 = (floor(i0), floor(j0), floor(k0))
            i1, j1, k1 = (i0 + 1, j0 + 1, k0 + 1)

            # Generate the grid point indices
            grid_indices = [(i, j, k) for k in (k0, k1) for j in (j0, j1)
                            for i in (i0, i1)]

            # Calculate coordinate at each of the grid index and calculate
            # the probability for each point
            calc_prob = lambda gp: self._calculate_point_weight(point, gp)
            grid_points = map(self.calculate_grid_coordinate, grid_indices)
            probability = map(calc_prob, grid_points)

            # Iterate through the grid indices and probability values and add
            # the count to the grid points around the given point.
            for (i, j, k), p in zip(grid_indices, probability):
                self.grid_data[grid_number, k, j, i] += count * p
