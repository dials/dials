from __future__ import absolute_import, division

def gz_open(filename, mode):
  import gzip
  return gzip.GzipFile(filename, mode)

def split_counts(image, split):
  from scitbx.array_family import flex
  import random

  new_images = [flex.int(flex.grid(image.focus()), 0) for k in range(split)]

  negative = image.as_1d() < 0
  positive = image.as_1d() > 0

  for new_image in new_images:
    new_image.as_1d().set_selected(negative, image.as_1d())

  for p in positive.iselection():
    counts = image[p]
    for j in range(counts):
      new_images[random.randint(0, split-1)][p] += 1

  return new_images

def merge_counts(images):
  from scitbx.array_family import flex
  image = flex.int(flex.grid(images[0].focus()), 0)
  negative = images[0].as_1d() < 0
  for i in images:
    image += i
  image.as_1d().set_selected(negative, images[0].as_1d())
  return image

def read_image(in_image):
  from scitbx.array_family import flex
  import binascii
  import os
  from dxtbx import load

  assert(os.path.exists(in_image))

  start_tag = binascii.unhexlify('0c1a04d5')

  data = gz_open(in_image, 'rb').read()
  data_offset = data.find(start_tag)
  cbf_header = data[:data_offset]
  pixel_values = load(in_image).get_raw_data()

  return pixel_values, cbf_header

def write_image(out_image, pixel_values, header):
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

  gz_open(out_image, 'wb').write(fixed_header + start_tag + compressed + tailer)

def main(in_images, out_images):
  assert(len(in_images) == len(out_images))
  n = len(in_images)
  import os
  for i in in_images:
    assert(os.path.exists(i))
  for o in out_images:
    assert(not os.path.exists(o))

  in_image_data = []
  in_image_headers = []

  for i in in_images:
    print "Reading %s" % i
    pixel, header = read_image(i)
    in_image_data.append(pixel)
    in_image_headers.append(header)

  sum_image = merge_counts(in_image_data)
  rebin_images = split_counts(sum_image, n)

  for o, pixel, header in zip(out_images, rebin_images, in_image_headers):
    print "Writing %s" % o
    write_image(o, pixel, header)
