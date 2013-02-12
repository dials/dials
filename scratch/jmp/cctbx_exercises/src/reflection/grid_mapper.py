
class GridMapper:
    """A class to map reflections from the detector to their own specific
    coordinate system grid.

    """

    def __init__(self, gonio, image, dcs, grid, wavelength, dp=(4, 4, 4)):
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
        self.image = image
        self.dcs = dcs
        self.grid = grid
        self.wavelength = wavelength
        self.dp = dp

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
                    s_dash = self.dcs.to_laboratory_coordinate(x, y)
                    s_dash = (1.0 / self.wavelength) * s_dash.normalize()
                    phi_dash = self.gonio.get_angle_from_frame(z)

                    # Calculate the reflection coordinate from the beam vector
                    c1, c2, c3 = rcs.from_diffracted_beam_vector(s_dash,
                                                                 phi_dash)

                    # Add the count from the pixel to the reflection grid
                    value = self.image[z-z0, y, x]
                    self.grid.add_point_count(index, (c1, c2, c3), value)
