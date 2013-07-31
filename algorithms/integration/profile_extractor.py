#
# copy_reflection_profiles.py
#
#  Copyright (C) 2013 Diamond Light Source
#
#  Author: James Parkhurst
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
from __future__ import division

class ProfileExtractor(object):
    ''' A class to extract the profiles from the sweep '''

    def __init__(self, sweep, crystal, mask=None, gain_map=None, dark_map=None,
                 **kwargs):
        ''' Initialise the class with the sweep etc

        Params:
            sweep The sweep to process
            crystal The crystal to process
            mask The detector mask
            gain_map The gain map
            dark_map The dark map
            **kwargs Other keyword arguments

        '''
        from dials.algorithms import shoebox
        from scitbx.array_family import flex

        # Save the sweep
        self.sweep = sweep

        # Set the detector mask
        image_size = sweep.get_image_size()[::-1]
        if not mask:
            self.mask = sweep[0] >= 0
        else:
            self.mask = mask

        # Set the gain map
        if not gain_map:
            self.gain_map = flex.double(flex.grid(*image_size), 1)
        else:
            self.gain_map = gain_map

        # Set the dark map
        if not dark_map:
            self.dark_map = flex.double(flex.grid(*image_size), 0)
        else:
            self.dark_map = dark_map

        # Get the parameters
        n_sigma = kwargs['bbox_nsigma']
        delta_d = n_sigma * sweep.get_beam().get_sigma_divergence(deg=False)
        delta_m = n_sigma * crystal.get_mosaicity(deg=False)

        # Create the function to mask the shoebox profiles
        self.mask_profiles = shoebox.Masker(sweep, mask, delta_d, delta_m)

    def __call__(self, reflections, adjacency_list=None):
        ''' Extract the profiles from the sweep

        Params:
            reflections The reflections to extract
            adjacency_list The list of overlapping reflections

        Returns:
            The reflection list

        '''
        from dials.algorithms import shoebox
        from dials.util.command_line import Command

        # Allocate memory for reflection profiles
        Command.start("Allocating reflection profiles")
        shoebox.allocate(reflections)
        Command.end("Allocated {0} reflection profiles".format(
            len([r for r in reflections if r.is_valid()])))

        # Mask the shoebox profiles
        self.mask_profiles(reflections, adjacency_list)

        # Extract the data from the sweep
        return self.extract_sweep_data(reflections)

    def extract_sweep_data(self, reflections):
        ''' Extract the profiles from the sweep

        Params:
            reflections The reflections to extract

        Returns:
            The reflection list

        '''
        from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
        from dials.util.command_line import ProgressBar
        from scitbx.array_family import flex
        from dials.algorithms import shoebox

        # Create the class to set all the shoebox pixels
        populate = shoebox.Populator(reflections, self.mask,
            self.gain_map, self.dark_map)

        # Create a progress bar
        progress = ProgressBar(title = "Extracting reflections")

        # For each image in the sweep, get the reflections predicted to have
        # been recorded on the image and copy the pixels from the image to
        # the reflection profile image. Update a progress bar as we go along
        first_array_index = self.sweep.get_array_range()[0]
        for index, image in enumerate(self.sweep):

            # Copy the image pixels from the image to the shoeboxes
            populate.add_image(image, index + first_array_index)
    #
    #        # Extract a mask from the shoeboxes showing which pixels to use
    #        mask = construct_image_mask_from_shoeboxes(detector_mask != 0,
    #          index + first_array_index, reflection_indices,
    #          reflections, kernel_size)
    #
    #        # Construct the threshold strategy
    #        threshold = XDSThresholdStrategy(
    #            kernel_size = kernel_size,
    #            gain = gain_map,
    #            mask = mask,
    #            n_sigma_b = n_sigma_b,
    #            n_sigma_s = n_sigma_s)
    #
    #        # Threshold the image
    #        mask = threshold(image)
    #
    #        # Assign the strong pixels and spots
    #        assign_strong_spots(mask, index + first_array_index,
    #            reflection_indices, reflections)
    #
            # Update the progress bar
            progress.update(100 * (index + 1) / len(self.sweep))

        # Progress bar finished
        progress.finished(
            "Extracted {0} reflections, {1} of which are strong".format(
                len([r for r in reflections if r.is_valid()]),
                len([r for r in reflections if r.is_strong()])))

        # Return the reflections
        return reflections
