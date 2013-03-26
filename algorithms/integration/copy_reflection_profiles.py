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
        f0 = r.bounding_box[0]
        f1 = r.bounding_box[1]
        for f in range(f0, f1):
            frames_to_reflection[f].append(i)

    # Return the list of lists
    return frames_to_reflection

def copy_image_pixels(sweep, reflections, frame_indices):
    """ Copy the image pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The list of reflections
        frame_indices The list of reflections recorded on each frame

    Returns:
        The updated list of reflections

    """
    from dials.algorithms.integration import copy_single_image_pixels
    from dials.util.command_line import ProgressBar

    # Create a progress bar
    progress = ProgressBar()

    # For each image in the sweep, get the reflections predicted to have
    # been recorded on the image and copy the pixels from the image to
    # the reflection profile image. Update a progress bar as we go along
    first_array_index = sweep.get_array_range()[0]
    for index, image in enumerate(sweep):
        reflection_indices = frame_indices[index]
        copy_single_image_pixels(image, index + first_array_index,
            reflection_indices, reflections)
        progress.update(100 * (index + 1) / len(sweep))

    # Progress bar finished
    progress.finished()

    # Return the reflections
    return reflections

def extract_reflection_profiles(sweep, reflections):
    """ Copy all the pixels from the sweep to the reflection profiles.

    Params:
        sweep The sweep object
        reflections The reflection list

    Returns:
        The updated reflection list.

    """
    from dials.algorithms.integration import allocate_reflection_profiles

    # Allocate memory for reflection profiles
    allocate_reflection_profiles(reflections)

    # Get the indices of the reflections recorded on each frame
    frame_indices = get_reflection_frame_indices(sweep, reflections)

    # Return the reflections
    return copy_image_pixels(sweep, reflections, frame_indices)
