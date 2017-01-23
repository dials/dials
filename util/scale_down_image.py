from __future__ import absolute_import, division

from dials.util import scale_down_array

def scale_down_array_py(image, scale_factor):
  '''Scale the data in image in a manner which retains the statistical structure
  of the input data. Input data type must be integers; negative values assumed
  to be flags of some kind (i.e. similar to Pilatus data) and hence preserved
  as input.'''

  from scitbx.random import variate, uniform_distribution
  from scitbx.array_family import flex

  assert (scale_factor <= 1)
  assert (scale_factor >= 0)

  # construct a random number generator in range [0, 1]
  dist = variate(uniform_distribution(0.0, 1.0))
  scaled_image = flex.int(len(image), 0)

  for j, pixel in enumerate(image):
    if pixel < 0:
      scaled_image[j] = pixel
    else:
      for c in range(pixel):
        if dist.next() < scale_factor:
          scaled_image[j] += 1

  return scaled_image

def read_image_to_flex_array(in_image):
  '''Looks like this works *only* for CBF images from a Pilatus detector;
  oh well - should still do something useful.'''
  from scitbx.array_family import flex
  import binascii
  import os
  from dxtbx import load

  assert(os.path.exists(in_image))

  start_tag = binascii.unhexlify('0c1a04d5')

  data = open(in_image, 'rb').read()
  data_offset = data.find(start_tag)
  cbf_header = data[:data_offset]
  pixel_values = load(in_image).get_raw_data()

  return pixel_values, cbf_header

def write_image_from_flex_array(out_image, pixel_values, header):
  '''Write a scaled CBF image from an array of pixel values and a header to
  add at the top. N.B. clobbers the binary size of the compressed data &
  the MD5 hash of the data.'''
  from cbflib_adaptbx import compress
  import binascii
  import os

  assert(not os.path.exists(out_image))
  start_tag = binascii.unhexlify('0c1a04d5')

  compressed = compress(pixel_values)

  fixed_header = ''
  for record in header.split('\n')[:-1]:
    if 'X-Binary-Size:' in record:
      old_size = int(record.split()[-1])
      fixed_header += 'X-Binary-Size: %d\r\n' % len(compressed)
    elif 'Content-MD5' in record:
      pass
    else:
      fixed_header += '%s\n' % record

  tailer = '\r\n--CIF-BINARY-FORMAT-SECTION----\r\n;\r\n'

  open(out_image, 'wb').write(fixed_header + start_tag + compressed + tailer)

def scale_down_image(in_image, out_image, scale_factor):
  '''Read in the data from in_image, apply the statistically valid scale factor
  to the data & write this out as out_image; retain the header as we go.'''

  image, header = read_image_to_flex_array(in_image)
  scaled_image = scale_down_array(image.as_1d(), scale_factor)
  from scitbx.array_family import flex
  sum_image = flex.sum(image.as_1d().select(image.as_1d() > 0))
  sum_scaled_image = flex.sum(scaled_image.select(scaled_image > 0))
  write_image_from_flex_array(out_image, scaled_image, header)
