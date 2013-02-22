

class Reflection:

    def __init__(self, hkl=None):

        self.hkl = hkl

    def calculate_phi_angles(self):
        pass

    def calculate_beam_vector(self):
        pass

    def calculate_detector_coordinates(self):
        pass

    def calculate_profile_coordinates(self):
        pass


if __name__ == '__main__':

    r = Reflection((0, 0, 0))
