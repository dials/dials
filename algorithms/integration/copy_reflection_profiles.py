from __future__ import division

def get_sweep_image_size(sweep):
    """Get the size of the image volume from the sweep."""

    # Get the detector image size
    detector_size = sweep.get_detector().get_image_size()

    # Get the scan range
    scan_range = sweep.get_scan().get_array_range()
    scan_size = scan_range[1] - scan_range[0]

    # return the image size
    return (scan_size, detector_size[1], detector_size[0])

def is_shoebox_valid(sweep, shoebox):
    """ Check if the shoebox is containined in the sweep."""

    # Get the image size
    image_size = get_sweep_image_size(sweep)

    # Ensure the shoebox is valid
    return (shoebox[0] >= 0 and shoebox[1] < image_size[0] and
            shoebox[2] >= 0 and shoebox[3] < image_size[1] and
            shoebox[2] >= 0 and shoebox[4] < image_size[2])

def extract_shoebox_pixels(sweep, shoebox):
    """Extract the image data from the sweep."""

    from scitbx.array_family import flex

    if is_shoebox_valid(sweep, shoebox):

      # Extract extents from shoebox
      si0, si1 = shoebox[0], shoebox[1]
      sj0, sj1 = shoebox[2], shoebox[3]
      sk0, sk1 = shoebox[4], shoebox[5]

      # Create the image array
      image = flex.int(flex.grid(sk1 - sk0, sj1 - sj0, si1 - si0))

      # Get the image size from the sweep
      image_size = get_sweep_image_size(sweep)

      # Copy all the pixels from the sweep to the image
      for k, sk in enumerate(range(sk0, sk1)):
          pixels = sweep[sk]
          for j, sj in enumerate(range(sj0, sj1)):
              for i, si in enumerate(range(si0, si1)):
                  image[k, j, i] = pixels[sj * image_size[0] + si]

    else:
        image = flex.int()

    # Return the image
    return image

def copy_reflection_profile(sweep, reflection):
    """Copy the reflection profile from the sweep to the reflection."""
    from sctibx.array_family import flex

    # Get the reflection image from the sweep
    reflection.image = extract_shoebox_pixels(sweep, reflection.shoebox)

    # Return the reflection
    return reflection

def copy_reflection_profiles(sweep, reflections):
    """Copy all the reflection profiles."""

    # Copy the profile image pixels for each reflection
    for r in reflections:
        copy_reflection_profile(r, sweep)

    # Return the reflections
    return reflections
