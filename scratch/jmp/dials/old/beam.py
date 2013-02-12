
class Beam:
    """Class to hold beam parameters"""

    def __init__(self, s0, wavelength):
        """Initialise the beam parameters

        Args:
            s0: The incident beam vector
            wavelength: The wavelength of the radiation

        """
        self.wavelength = wavelength
        self.s0 = s0.normalize() / self.wavelength
