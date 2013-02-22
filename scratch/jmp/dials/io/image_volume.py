def get_image_volume_size(list_of_images):
    """Get the size of the image volume"""
    import os
    import numpy
    from dxtbx.model.detector_helpers_types import detector_helpers_types
    from dxtbx.format.Registry import Registry

    for image in list_of_images:
        assert(os.path.exists(image))

    list_of_images.sort()

    format = Registry.find(list_of_images[0])

    # verify that these are all the same format i.e. that they are all
    # understood equally well by the format instance.

    format_score = format.understand(list_of_images[0])

    for image in list_of_images:
        assert(format.understand(image) == format_score)

    i = format(list_of_images[0])

    beam = i.get_beam()
    gonio = i.get_goniometer()
    det = i.get_detector()
    scan = i.get_scan()

    # now verify that they share the same detector position, rotation axis
    # and beam properties.

    scans = [scan]

    for image in list_of_images[1:]:
        i = format(image)
        assert(beam == i.get_beam())
        assert(gonio == i.get_goniometer())
        assert(det == i.get_detector())

        scans.append(i.get_scan())

    for s in sorted(scans)[1:]:
        scan += s

    size_z = scan.get_image_range()[1] - scan.get_image_range()[0]
    size_xy = det.get_image_size()

    return (size_z, size_xy[1], size_xy[0])

def get_image_volume(image_paths):
    """Get the image volume from the given list of files"""
    from iotbx.detectors import ImageFactory
    import numpy
    num_slices, height, width = get_image_volume_size(image_paths)

    # Initialise the image volume
    num_slices = len(image_paths)
    volume = numpy.zeros(shape=(num_slices, height, width), dtype=numpy.int32)

    # For each CBF file, read the image and put into the image volume
    for i, filename in enumerate(image_paths):
        image = ImageFactory(filename)
        image.read()
        image_data = image.get_raw_data().as_numpy_array()
        image_data.shape = (height, width)
        volume[i,:,:] = image_data

    # Return image volume
    return volume

def search_for_image_volume(search_path):
    """Search for the image volume"""
    from glob import glob
    filenames = glob(search_path)
    filenames.sort()
    return get_image_volume(filenames)
