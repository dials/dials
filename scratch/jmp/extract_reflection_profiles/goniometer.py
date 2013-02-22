
class CoordinateSystem:
    """The goniometer coordinate system"""

    def __init__(self, m2, s0):
        """Create the goniometer coordinate system

        Args:
            m2: The rotation axis
            s0: The incident beam vector

        """

        # Set the incident beam vector and rotation axis
        self.s0 = s0
        self.m2 = m2

        # Calculate the other goniometer axes
        self.m1 = self.m2.cross(self.s0).normalize()
        self.m3 = self.m1.cross(self.m2).normalize()

    def get_rotation_axis(self):
        """Get the rotation axis

        Returns:
            The rotation axis

        """
        return self.m2


class Goniometer:
    """The goniometer parameters"""

    def __init__(self, m2, s0, z0, phi0, dphi):
        """Initialise the goniometer parameters.

        Args:
            z0: The starting frame
            phi0: The starting angle
            dphi: The oscillation range

        """
        # Create the goniometer coordinate system
        self.coordinate_system = CoordinateSystem(m2, s0)

        # Set the other goniometer parameters
        self.z0 = z0
        self.phi0 = phi0
        self.dphi = dphi

    def get_angle_from_frame(self, z):
        """Calculate the angle from the frame number.

        Args:
            z: The frame number

        Returns:
            The phi angle at that frame number

        """
        return self.phi0 + (z - self.z0) * self.dphi

    def get_frame_from_angle(self, phi):
        """Calculate the frame number at the given angle.

        The number returned is a float.

        Args:
            phi: The phi angle

        Returns:
            The frame number at that angle.

        """
        return self.z0 + (self.phi - self.phi0) / self.dphi
