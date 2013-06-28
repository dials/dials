from __future__ import division
class ReflectionMaskCreator(object):
    """Class to create the reflection mask."""

    def __init__(self, beam, detector, gonio, detector_mask, volume_size,
                 sigma_divergence, sigma_mosaicity, n_sigma = 10,
                 roi_volume_fraction = 0.99):
        """Setup the reflection mask creator.

        Args:
            beam The beam parameters
            detector The detector parameters
            gonio The goniometer parameters
            detector_mask The mask of bad pixels on the detector
            volume_size The size of the image volume
            sigma_divergence The standard deviation of the beam divergence
            sigma_mosaicity The standard deviation of the mosaicity
            n_sigma The number of standard deviations / 2 to use

        """
        from dials_jmp.integration import ReflectionMaskRoi, ReflectionMask

        # Create the reflection mask roi calculator
        self.reflection_mask_roi = ReflectionMaskRoi(
                                    beam,
                                    detector,
                                    gonio,
                                    n_sigma * sigma_divergence,
                                    n_sigma * sigma_mosaicity)

        # Create the reflection mask
        self.reflection_mask = ReflectionMask(detector_mask, volume_size)

        # Save the roi volume fraction
        self.roi_volume_fraction = roi_volume_fraction

    def create(self, reflections):
        """Create the reflection mask.

        First calculate the roi for each reflections. Then remove reflections
        with an roi volume above a certain percentage of the maximum volume.
        Then create the reflection mask volume and filter out any further
        reflections that are invalid.

        Args:
            reflections The list of reflections

        Returns:
            The modified list of reflections

        """
        from dials_jmp.array_family import flex
        from dials_jmp.integration import filter_reflections_by_roi_volume

        # Create the reflection region of interests
        self.reflection_mask_roi.calculate(reflections)

        # Filter the reflections by region of interest volume
        reflections = flex.remove_if_not(reflections,
                        filter_reflections_by_roi_volume(
                            reflections, self.roi_volume_fraction))

        # Create the reflection mask itself and filter any bad reflections
        reflections = flex.remove_if_not(reflections,
                        self.reflection_mask.create(reflections))

        # Return the reflections
        return reflections

    def set_reflection_pixel_value(self, image_volume, reflection, value):
        """Set the image volume value for the roi"""
        roi = reflection.region_of_interest
        for k in range(roi[4], roi[5]):
            for j in range(roi[2], roi[3]):
                for i in range(roi[0], roi[1]):
                    if (self.mask[k, j, i] == reflection.mask_index):
                        image_volume[k, j, i] = value

    @property
    def mask(self):
        """Return the reflection mask array"""
        return self.reflection_mask.mask
