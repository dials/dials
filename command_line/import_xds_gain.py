from __future__ import division

def import_xds_gain(sweep_filename, in_filename, out_filename):
  from dials.util import image
  from dials.model.serialize import load
  from math import ceil

  # Read the sweep
  print "Reading sweep from %s" % sweep_filename
  sweep = load.sweep(sweep_filename)

  # Read the image in
  print "Reading gain image from %s" % in_filename
  inhandle = image.reader()
  inhandle.read_file(in_filename)
  sampled_gain = inhandle.get_data().as_numpy_array()

  # Create the output image
  assert(len(sweep.get_detector()) == 1)
  image_size = sweep.get_detector()[0].get_image_size()[::-1]

  # Resize the gain image
  print "Creating new image"
  import numpy
  ny = int(ceil(image_size[1] / sampled_gain.shape[0]))
  nx = int(ceil(image_size[0] / sampled_gain.shape[1]))
  i1 = numpy.array([range(image_size[1])] * image_size[0], dtype=numpy.int32)
  j1 = numpy.array([range(image_size[0])] * image_size[1], dtype=numpy.int32).transpose()
  i2 = numpy.divide(i1, nx)
  j2 = numpy.divide(j1, ny)
  gain = numpy.zeros(image_size, dtype=numpy.float64)
  gain[j1,i1] = sampled_gain[j2,i2] / 1000.0

  # Write the file to CBF
  print "Writing gain image to %s" % out_filename

if __name__ == '__main__':
  from optparse import OptionParser
  usage = "usage: %prog [options] /path/to/sweep.json /path/to/gain/file"
  parser = OptionParser(usage)

  # Write the datablock to JSON or Pickle
  parser.add_option(
    "-o", "--output",
    dest = "output",
    type = "string", default = "dials_gain.cbf",
    help = "The output gain file")

  # Parse the command line arguments
  (options, args) = parser.parse_args()
  if len(args) < 2:
    parser.print_help()

  # Import the gain file
  import_xds_gain(args[0], args[1], options.output)
