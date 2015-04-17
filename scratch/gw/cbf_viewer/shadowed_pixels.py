# -*- coding: utf-8 -*-

"""
Contains the Pixshadow class, used to calculate shadow mask.
"""

import analyser
import shadowmapper

class Pixshadow(object):
    """
    Uses shadowmapper to return indices and affected pixels.

    args:
    filename: The cbf image to base the detector on. The detector
    dimensions and position will be extracted from this.
    """
    def __init__(self, filename):
        self.newimage = analyser.Image(filename)

        (origin, fast, slow, dimensions, pix, dist,
        res) = self.newimage.detector_parameters()

        self.newshadow = shadowmapper.Shadow(shadowmapper.day2_coords,
                          filename, det_angle, dist, origin, fast, slow)

    def pixcounter(self, kappa, omega):
        """
        Compares a shadowed image to an unshadowed image.

        Uses the shadowmapper.Shadow class's affected_pixels method,
        returning the number of affected pixels and their indices.

        args:
        kappa, omega: The goniometer angles to use.
        """
        affected, indices = self.newshadow.affected_pixels(kappa, omega,
                                                              True)[1:3]
        return affected, indices


if __name__ == "__main__":
    pix = Pixshadow('cbf_image2.cbf')
    a, b = pix.pixcounter(-180.0, 90.0)
    print a
    print b
