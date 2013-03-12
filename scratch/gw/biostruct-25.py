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

if __name__ == '__main__':
    import sys
    image = open_file_return_array(sys.argv[1])
    plot_image(image)
