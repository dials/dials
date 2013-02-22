
class CoordinateSystem:
    """ A class to represent the detector coordinate system. The coordinate
    system is defined, as in XDS, in such a way that the laboratory coordinate
    at a given pixel, (x, y) is as follows:

        (x - x0)d1 + (y - y0)d2 + fd3

    Where:
        x0, y0 is the origin of the detector in pixel coordinates
        d1, d2 are the x and y detector axis vectors
        d3 is the detector normal
        f is the distance to the detector

    """

    def __init__(self, d1, d2, d3, distance, pixel_origin, pixel_size):
        """Initialise variables for coordinate transform between the laboratory
        and detector coordinate systems.

        Args:
            d1: The detector x axis unit vector
            d2: The detector y axis unit vector
            d3: The detector normal unit vector
            distance: The distance from the crystal to the detector
            pixel_origin: The (x0, y0) pixel coordinates of the detector origin
            pixel_size: The size of the pixels

        """
        # Save the detector distance, the pixel sizes and the pixel origin
        self.f = distance
        self.px, self.py = pixel_size
        self.x0, self.y0 = pixel_origin

        # Save the axis unit vectors
        self.d1 = d1.normalize()
        self.d2 = d2.normalize()
        self.d3 = d3.normalize()

        # The d1, d2 and d3 given should be unit vectors. In order to transform
        # between laboratory and detector (pixel) coordinates, we need to have
        # two sets of axes vectors. One set scaled by the pixel size and the
        # other scales by 1 / (pixel size).
        self.d1_0 = self.d1 * self.px
        self.d2_0 = self.d2 * self.py
        self.d1_1 = self.d1 / self.px
        self.d2_1 = self.d2 / self.py

        # Calculate the origin in laboratory coordinates
        self.origin = self.f * self.d3

    def get_axes(self):
        """Get a tuple containing the detector axes

        Returns:
            A tuple containing the detector axes.

        """
        return tuple((self.d1, self.d2, self.d3))

    def to_laboratory_coordinate(self, x, y):
        """Calculate the laboratory coordinates of a detector pixel

        Args:
            x: The x detector pixel coordinate
            y: The y detector pixel coordinate

        Returns:
            The laboratory coordinate

        """
        return (x - self.x0) * self.d1_0 + \
               (y - self.y0) * self.d2_0 + self.origin

    def from_diffracted_beam_vector(self, s):
        """Calculate the detector pixel coordinate from the beam vector

        Args:
            s: The beam vector in laboratory coordinates

        Returns:
            The (x, y) pixel coordinate or None.

        """
        # Calculate the point (in pixels) at which the diffracted beam vector
        # intersects the detector plane. (only when f * s1.d3 > 0)
        s_dot_d3 = s.dot(self.d3)
        if self.f * s_dot_d3 <= 0:
            return (None, None)

        x = self.x0 + self.f * s.dot(self.d1_1) / s_dot_d3
        y = self.y0 + self.f * s.dot(self.d2_1) / s_dot_d3
        return (x, y)
