from __future__ import division

class CoordinateSystem:
    """The profile specific coordinate system."""

    def __init__(self, s0, s1, m2, phi):
        """Construct the profile coordinate system.

        Args:
            s0: The incident beam vector
            s1: The diffracted beam vector
            m2: The rotation axis
            phi: The rotation angle

        """
        from math import pi

        # Save the input variables
        self.s0 = s0
        self.s1 = s1
        self.m2 = m2
        self.phi = phi

        # Calculate p*
        self.p_star = self.s1 - self.s0

        # Calculate the scale factors
        self.scale  = 180.0 / (pi * self.s1.length())
        self.scale3 = 180.0 / (pi * self.p_star.length())

        # Calculate the three axes
        self.e1 = self.s1.cross(self.s0).normalize()
        self.e2 = self.s1.cross(self.e1).normalize()
        self.e3 = (self.s1 + self.s0).normalize()

        # Calculate zeta  (related to the lorentz correction factor)
        self.zeta = m2.dot(self.e1)

    def get_axes(self):
        """Get the coordinate system axes as a tuple.

        Returns:
            The coordinate system as a tuple (e1, e2, e3)

        """
        return (tuple(self.e1), tuple(self.e2), tuple(self.e3))

    def from_diffracted_beam_vector(self, s_dash, phi_dash):
        """Calculate a coordinate in the profile system

        Args:
            s_dash: A diffracted beam vector
            phi_dash: The beam vector's rotation angle

        Returns:
            A tuple, (x1, x2, x3), containing the point's coordinates.

        """
        # Calculate the rotation of p_star about m2 with angle of dphi
        p_star_r = self.p_star.rotate(self.m2, phi_dash - self.phi, True)

        # Calculate the three points
        x1 = self.e1.dot(s_dash - self.s1) * self.scale
        x2 = self.e2.dot(s_dash - self.s1) * self.scale
        x3 = self.e3.dot(p_star_r - self.p_star) * self.scale3
        #x3 = self.zeta * (phi_dash - self.phi)
        return (x1, x2, x3)


if __name__ == '__main__':

    from scitbx import matrix

    wavelength = 0.689400
    s0 = matrix.col(( 0.013142,  0.002200,  1.450476)).normalize() / wavelength
    s1 = matrix.col((-0.017520, -0.247753,  1.428447)).normalize() / wavelength
    m2 = matrix.col(( 0.999975, -0.001289, -0.006968))
    phi = 5.835757

    cs = CoordinateSystem(s0, s1, m2, phi)

    e1 = cs.e1
    e2 = cs.e2
    e3 = cs.e3

    phi_dash = int(phi)
    s_dash = s1.rotate(m2, phi - phi_dash, True)
    #s_dash = s1

    print cs.calculate_coordinate(s_dash, phi_dash)
