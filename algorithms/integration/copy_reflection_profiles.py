from __future__ import division


def get_reflections_recorded_on_frames(sweep, reflections):

    from collections import defaultdict

    # Create a list of lists
    frames_to_reflection = defaultdict(list)

    # For each reflection, Find the frames which it spans and copy an
    # index into the frame -> reflection list
    for i, r in enumerate(reflections):
        f0 = r.shoebox[4]
        f1 = r.shoebox[5]
        for f in range(f0, f1):
            frames_to_reflection[f].append(i)

    # Return the list of lists
    return frames_to_reflection


def is_shoebox_valid(sweep, shoebox):
    """ Check if the shoebox is containined in the sweep."""

    size = (len(sweep), sweep.get_image_size()[1], sweep.get_image_size()[0])

    # Ensure the shoebox is valid
    return (shoebox[0] >= 0 and shoebox[1] < size[2] and
            shoebox[2] >= 0 and shoebox[3] < size[1] and
            shoebox[4] >= 0 and shoebox[5] < size[0])

def extract_shoebox_pixels(sweep, shoebox):
    """Extract the image data from the sweep."""

    from scitbx.array_family import flex

    if is_shoebox_valid(sweep, shoebox):

      # Extract extents from shoebox
      si0, si1 = shoebox[0], shoebox[1]
      sj0, sj1 = shoebox[2], shoebox[3]
      sk0, sk1 = shoebox[4], shoebox[5]

      # Get the image data from the sweep
      image = sweep.to_array((sk0, sk1, sj0, sj1, si0, si1))

    else:
        image = flex.int()

    # Return the image
    return image

def copy_reflection_profile(sweep, reflection):
    """Copy the reflection profile from the sweep to the reflection."""
    from scitbx.array_family import flex

    # Get the reflection image from the sweep
    reflection.image = extract_shoebox_pixels(sweep, reflection.shoebox)

    # Return the reflection
    return reflection

def copy_reflection_profiles(sweep, reflections):
    """Copy all the reflection profiles."""
    import numpy
    from scitbx.array_family import flex

    # Get an index array for the reflection list
    indices = range(len(reflections))

    # Sort the images by z
    print "Sorting"
    indices = sorted(indices, key=lambda x: reflections[x].frame_number)

    print "Getting frames"
    frames_to_reflections = get_reflections_recorded_on_frames(
        sweep, reflections)


    # Allocate all images
    print "Allocating"
    for i in range(len(reflections)):
        r = reflections[i]
        size_x = r.shoebox[1] - r.shoebox[0]
        size_y = r.shoebox[3] - r.shoebox[2]
        size_z = r.shoebox[5] - r.shoebox[4]
        r.image = flex.int(flex.grid(size_z, size_y, size_x))
        reflections[i] = r


    image_size = sweep.get_image_size()


    for index, image in enumerate(sweep):
        print index
        rr = frames_to_reflections[index]
        for ri in rr:

            r = reflections[ri]

            # Extract extents from shoebox
            si0, si1 = r.shoebox[0], r.shoebox[1]
            sj0, sj1 = r.shoebox[2], r.shoebox[3]
            sk0, sk1 = r.shoebox[4], r.shoebox[5]
            sk = index - sk0
            if si0 >= 0 and sj0 >= 0 and si1 < image_size[0] and sj1 < image_size[1]:
                npa = r.image.as_numpy_array()

                npa[sk,:,:] = image[sj0:sj1,si0:si1].as_numpy_array()
                r.image = flex.int(npa)
                reflections[ri] = r

    # Return the reflections
    return reflections
