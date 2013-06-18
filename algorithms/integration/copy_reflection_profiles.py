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

def get_reflection_frame_indices(sweep, reflections):
    """ Go through the reflection list and check for each frame which
    reflections are recorded on which frames as defined by the reflection
    shoebox.

    Params:
        sweep The sweep object
        reflections The reflection list

    Returns:
        A list of frames with a list of recorded reflections

    """
    from collections import defaultdict
    from scitbx.array_family import flex

    # Create a dictionary of flex arrays
    frames_to_reflection = defaultdict(flex.int)

    # For each reflection, Find the frames which it spans and copy an
    # index into the frame -> reflection list
    for i, r in enumerate(reflections):
        f0 = r.bounding_box[4]
        f1 = r.bounding_box[5]
        for f in range(f0, f1):
            frames_to_reflection[f].append(i)

    # Return the list of lists
    return frames_to_reflection

def copy_image_pixels(sweep, reflections, frame_indices,
                      gain_map=None, dark_map=None):
    """ Copy the image pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The list of reflections
        frame_indices The list of reflections recorded on each frame
        gain_map The gain map
        dark_map The dark map

    Returns:
        The updated list of reflections

    """
    from dials.algorithms.integration import copy_single_image_pixels
    from dials.util.command_line import ProgressBar
    from scitbx.array_family import flex

    # If gain map or dark map are None then set suitable defaults
    image_size = sweep.get_image_size()[::-1]
    if not gain_map:
        gain_map = flex.double(flex.grid(*image_size), 1)
    if not dark_map:
        dark_map = flex.int(flex.grid(*image_size), 0)

    # Create a progress bar
    progress = ProgressBar(title="Extracting reflections")

    # For each image in the sweep, get the reflections predicted to have
    # been recorded on the image and copy the pixels from the image to
    # the reflection profile image. Update a progress bar as we go along
    first_array_index = sweep.get_array_range()[0]
    for index, image in enumerate(sweep):
        reflection_indices = frame_indices[index]
        copy_single_image_pixels(image, index + first_array_index,
            reflection_indices, reflections, gain_map, dark_map)
        progress.update(100 * (index + 1) / len(sweep))

    # Progress bar finished
    progress.finished("Extracted {0} reflections".format(len(reflections)))

    # Return the reflections
    return reflections

def extract_reflection_profiles(sweep, reflections, adjacency_list=None,
                                gain_map=None, dark_map=None):
    """ Copy all the pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The reflection list
        adjacency_list The adjacency list (optional)
        gain_map The detector gain map
        dark_map The detector dark map

    Returns:
        The updated reflection list.

    """
    from dials.algorithms.integration import allocate_reflection_profiles
    from dials.algorithms.integration import ShoeboxMasker
    from scitbx.array_family import flex
    from dials.util.command_line import Command

    # Allocate memory for reflection profiles
    Command.start("Allocating reflection profiles")
    reflections = allocate_reflection_profiles(reflections)
    Command.end("Allocated {0} reflection profiles".format(len(reflections)))

    # Get the indices of the reflections recorded on each frame
    Command.start("Getting reflection frame indices")
    frame_indices = get_reflection_frame_indices(sweep, reflections)
    Command.end("Got frame indices for {0} reflections".format(len(reflections)))

    # Copy the pixels from the sweep to the reflection shoeboxes
    reflections = copy_image_pixels(sweep, reflections, frame_indices,
        gain_map, dark_map)

    # If the adjacency list is given, then create the reflection mask
    if adjacency_list:
        detector_mask = (sweep[0] >= 0).as_1d().as_int()
        detector_mask.reshape(flex.grid(sweep[0].all()))
        Command.start("Masking overlapped reflections")
        shoebox_masker = ShoeboxMasker(detector_mask)
        shoebox_masker(reflections, adjacency_list)
        Command.end("Masked {0} overlapped reflections".format(
            len(adjacency_list)))

    # Return the reflections
    return reflections
