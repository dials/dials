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

def copy_image_pixels(sweep, reflections, gain_map = None, dark_map = None,
                      kernel_size = (3, 3), n_sigma_b = 6.0, n_sigma_s = 3.0,
                      detector_mask = None):
    """ Copy the image pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The list of reflections
        frame_indices The list of reflections recorded on each frame
        gain_map The gain map
        dark_map The dark map
        kernel_size The size of the smoothing kernel to use
        n_sigma_b The number of sigmas for the background
        n_sigma_s The number of sigmas for a strong pixel
        detector_mask The detector mask

    Returns:
        The updated list of reflections

    """
    from dials.algorithms.integration import assign_strong_spots
    from dials.algorithms.peak_finding.threshold import XDSThresholdStrategy
    from dials.util.command_line import ProgressBar, Command
    from scitbx.array_family import flex
    from dials.algorithms import shoebox

    # If gain map or dark map are None then set suitable defaults
    image_size = sweep.get_image_size()[::-1]
    if not gain_map:
        gain_map = flex.double(flex.grid(*image_size), 1)
    if not dark_map:
        dark_map = flex.double(flex.grid(*image_size), 0)

    # Set the detector mask if none is set
    if detector_mask == None:
        detector_mask = sweep[0] >= 0

    # Create the class to set all the shoebox pixels
    populate = shoebox.Populator(reflections, detector_mask, gain_map, dark_map)

    # Create a progress bar
    progress = ProgressBar(title = "Extracting reflections")

    # For each image in the sweep, get the reflections predicted to have
    # been recorded on the image and copy the pixels from the image to
    # the reflection profile image. Update a progress bar as we go along
    first_array_index = sweep.get_array_range()[0]
    for index, image in enumerate(sweep):

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
        progress.update(100 * (index + 1) / len(sweep))

    # Progress bar finished
    progress.finished(
        "Extracted {0} reflections, {1} of which are strong".format(
            len([r for r in reflections if r.is_valid()]),
            len([r for r in reflections if r.is_strong()])))

    # Return the reflections
    return reflections

def extract_reflection_profiles(sweep, reflections, adjacency_list = None,
                                gain_map = None, dark_map = None, kernel_size = (3, 3),
                                n_sigma_b = 6.0, n_sigma_s = 3.0,
                                detector_mask = None, delta_d=None, delta_m=None):
    """ Copy all the pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The reflection list
        adjacency_list The adjacency list (optional)
        gain_map The detector gain map
        dark_map The detector dark map
        kernel_size The size of the smoothing kernel to use
        n_sigma_b The number of sigmas for the background
        n_sigma_s The number of sigmas for a strong pixel
        detector_mask The detector mask

    Returns:
        The updated reflection list.

    """
    from dials.algorithms import shoebox
    from scitbx.array_family import flex
    from dials.util.command_line import Command

    # Allocate memory for reflection profiles
    Command.start("Allocating reflection profiles")
    shoebox.allocate(reflections)
    Command.end("Allocated {0} reflection profiles".format(
        len([r for r in reflections if r.is_valid()])))

    # If the adjacency list is given, then create the reflection mask
    detector_mask = (sweep[0] >= 0)
    shoebox_masker = shoebox.masker.Masker(sweep, detector_mask,
        delta_d, delta_m)
    shoebox_masker(reflections, adjacency_list)


    # Copy the pixels from the sweep to the reflection shoeboxes
    reflections = copy_image_pixels(sweep, reflections,
        gain_map, dark_map, kernel_size, n_sigma_b, n_sigma_s,
        detector_mask)

    # Return the reflections
    return reflections
