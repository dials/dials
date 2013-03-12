from __future__ import division
# copied from James's code...

def open_file_return_array(filename):
    import pycbf
    import numpy

    # open file
    cbf_handle = pycbf.cbf_handle_struct()
    cbf_handle.read_file(filename, pycbf.MSG_DIGEST)
    cbf_handle.rewind_datablock()

    # find right datablock
    cbf_handle.select_datablock(0)
    cbf_handle.select_category(0)
    cbf_handle.select_column(2)
    cbf_handle.select_row(0)

    type = cbf_handle.get_typeofvalue()
    assert (type.find('bnry') > -1)

    # read and reshape
    image = numpy.fromstring(cbf_handle.get_integerarray_as_string(),
                             numpy.int32)
    parameters = cbf_handle.get_integerarrayparameters_wdims()
    image.shape = (parameters[10], parameters[9])

    return image

def plot_image(image):
    from matplotlib import pyplot
    pyplot.imshow(image)
    pyplot.savefig('biostruct-25.png')

    return

def read_integrate_hkl_apply_corrections(int_hkl, x_crns_file, y_crns_file):
    '''Read X corrections and Y corrections to arrays; for record in
    int_hkl read x, y positions and compute r.m.s. deviations between
    observed and calculated centroids with and without the corrections.
    N.B. will only compute offsets for reflections with I/sigma > 10.'''

    from scitbx import matrix
    import math

    x_corrections = open_file_return_array(x_crns_file)
    y_corrections = open_file_return_array(y_crns_file)

    sumxx_orig = 0.0
    n_obs = 0
    sumxx_corr = 0.0

    for record in open(int_hkl):
        if record.startswith('!'):
            continue
        values = map(float, record.split())
        i, sigi = values[3], values[4]

        if sigi < 0:
            continue

        if i / sigi < 10:
            continue

        xc, yc, zc = values[5:8]
        xo, yo, zo = values[12:15]

        xc_orig = matrix.col((xc, yc))
        xo_orig = matrix.col((xo, yo))

        n_obs += 1
        sumxx_orig += (xo_orig - xc_orig).dot()

        # get the correcions for this pixel... N.B. 1 + ... lost due to C-style
        # rather than Fortran-style arrays

        ix4 = int(round((xc - 2) / 4))
        iy4 = int(round((yc - 2) / 4))

        # hmm.... Fortran multi-dimensional array ordering... identified by
        # bounds error other way around...

        dx = 0.1 * x_corrections[iy4, ix4]
        dy = 0.1 * y_corrections[iy4, ix4]

        xc_corr = matrix.col((xc + dx, yc + dy))

        sumxx_corr += (xo_orig - xc_corr).dot()

    return math.sqrt(sumxx_orig / n_obs), math.sqrt(sumxx_corr / n_obs)

if __name__ == '__main__':
    import sys

    orig, corr = read_integrate_hkl_apply_corrections(
        sys.argv[1], sys.argv[2], sys.argv[3])

    print orig, corr
