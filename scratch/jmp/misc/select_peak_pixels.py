from __future__ import division

def select_contigious(index, shape, conn):
    """ Filter the indices to only contain contigious pixels

    Params:
        index The input array of indices
        shape The shape of the mask
        conn The connectivity (4 or 8)

    Returns:
        The filtered indices

    """
    from scipy.ndimage.measurements import label, histogram
    from numpy import zeros, int32, argmax, where
    from operator import mul

    # Set the label structure
    if conn == 8:
        structure = [[1, 1, 1], [1, 1, 1], [1, 1, 1]]
    else:
        structure = [[0, 1, 0], [1, 1, 1], [0, 1, 0]]

    # Create a mask of the roi size
    mask = zeros(reduce(mul, shape), dtype=int32)

    # Set all the points at the indices to 1
    for i in index:
      mask[i] = 1

    # Reshape the mask
    mask.shape = shape

    # Label the indices in the mask
    regions, nregions = label(mask, structure)

    # Histogram the regions to get the largest
    histo = histogram(regions, 0, nregions+1, nregions+1)

    # Find largest region
    max_index = argmax(histo[1:]) + 1

    # Return only those region indices that are equal to max_index
    regions.shape = (-1)
    return flex.size_t(where(regions == max_index)[0])


def select_foreground_pixels(pixel_data, min_pixels=10, n_sigma=-1, conn=4):
    """ Select the pixels contributing to the spot foreground

    Fixed select all foreground pixels that leave the remaining background
    pixels in an approximate normal distribution. Then filter indices to
    get only those pixels that are contigious

    Params:
        pixel_data The pixel data
        min_pixels The minimum number of pixels contributing to background
        n_sigma The number of standard deviations to check for
        conn The connectivity (4 or 8)

    Returns:
        A list of indices deemed to contribute to spot foreground

    """
    from dials.algorithms.integration import foreground_pixels
    from scitbx.array_family import flex

    # Get the 1D array
    data = flex.double(len(pixel_data))
    for i in range(len(pixel_data)):
        data[i] = pixel_data[i]

    # Get a list of foreground pixel indices
    index = foreground_pixels(data, min_pixels, n_sigma)

    # Select and return only those pixels that are contigious
    return select_contigious(index, (5, 5), conn)




if __name__ == '__main__':

    from scitbx.array_family import flex
    data = flex.double([[1, 10,  0,  0, 0],
                        [0,  0,  0,  1, 1],
                        [2,  1, 10,  0, 0],
                        [0,  0, 10, 20, 0],
                        [0,  0,  1,  0, 1]])

    index = select_foreground_pixels(data, n_sigma=2)

    for i in index:
        print i
