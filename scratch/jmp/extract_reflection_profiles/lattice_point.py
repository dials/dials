from __future__ import division

class LatticePoint:
    """A class to represent a reciprocal lattice point"""

    def __init__(self, hkl):
        """Initialise the lattice point with the miller indices.

        Args:
            hkl: The (h, k, l) miller indices

        """
        self.hkl = hkl
        self.phi = ()

    def calculate_intersection_angles(self, ra):
        """From a previously specified rotation angle object, calculate the
        intersection angles. If sucessful, two will be generated and the
        internal self.phi array will be set and the function will return True,
        otherwise return False.

        Args:
            ra: The rotation angles object

        Returns:
            True/False

        """
        import math

        # Convert radians to degrees
        r2d = 180.0 / math.pi

        # Calculate the 2 intersection angles for the reflection index
        # and append them to the list of phi-angles
        if ra(self.hkl):
            angles = ra.get_intersection_angles()
            self.phi = (angles[0] * r2d % 360, angles[1] * r2d % 360)
            return True

        # Return failed
        return  False
