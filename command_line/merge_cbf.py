#!/usr/bin/env python
# merge_cbf.py
#
#   Copyright (C) 2015 Diamond Light Source, Richard Gildea
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#

from __future__ import division

import iotbx.phil

help_message = '''

Examples::

  dials.merge_cbf image_*.cbf 10

'''

phil_scope = iotbx.phil.parse("""\
merge_n_images = 2
  .type = int(value_min=1)
  .help = "Number of input images to average into a single output image"
output {
  image_prefix = sum_
    .type = path
}
""", process_includes=True)

def merge_cbf(imageset, n_images, out_prefix="sum_"):

  from dxtbx.format.FormatCBF import FormatCBF
  assert issubclass(imageset.reader().get_format_class(), FormatCBF), (
    "Only CBF format images supported")

  from cbflib_adaptbx import compress
  import binascii

  assert len(imageset) >= n_images

  n_output_images = len(imageset) // n_images

  in_oscillation = imageset.get_scan().get_oscillation()[1]
  out_oscillation = in_oscillation * n_images

  for i_out in range(n_output_images):
    data_out = None

    for j in range(n_images):

      i_in = (i_out*n_images) + j

      data_in = imageset[i_in]
      if data_out is None:
        data_out = data_in
      else:
        data_out += data_in

    out_image = "%s%04i.cbf" %(out_prefix, i_out+1)

    start_tag = binascii.unhexlify('0c1a04d5')

    data = open(imageset.get_path(i_out*n_images), 'rb').read()
    data_offset = data.find(start_tag)
    cbf_header = data[:data_offset]

    new_header = []
    compressed = compress(data_out)

    old_size = 0

    for record in cbf_header.split('\n')[:-1]:
      if 'X-Binary-Size:' in record:
        old_size = int(record.split()[-1])
        new_header.append('X-Binary-Size: %d\r\n' % len(compressed))
      elif 'Content-MD5' in record:
        pass
      elif '# Angle_increment' in record:
        new_header.append('# Angle_increment %.4f deg.\r\n' %out_oscillation)
      else:
        new_header.append('%s\n' % record)

    tailer = data[data_offset + 4 + old_size:]

    with open(out_image, 'wb') as f:
      f.write(''.join(new_header) + start_tag + compressed + tailer)
      print '%s written' % out_image

  return

if __name__ == '__main__':
  import sys
  import libtbx.load_env

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks

  usage = "%s [options] image_*.cbf" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    epilog=help
  )

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)

  n_images = params.merge_n_images
  out_prefix = params.output.image_prefix
  datablocks = flatten_datablocks(params.input.datablock)
  if len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]

  merge_cbf(imageset, n_images, out_prefix=out_prefix)
