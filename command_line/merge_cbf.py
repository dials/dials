#!/usr/bin/env python
# merge_cbf.py
#
#   Copyright (C) 2015 Diamond Light Source, Richard Gildea
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#

from __future__ import absolute_import, division

import iotbx.phil

help_message = '''

This program can be used to merge a given number of consecutive cbf files into
a smaller number of images For example, running dials.merge_cbf on a datablock
with 100 images, using the default value of merge_n_images=2, will output 50
summed images, with every consecutive pair of images being summed into a single
output image. Currently only cbf format images are supported as input.

Examples::

  dials.merge_cbf image_*.cbf

  dials.merge_cbf image_*.cbf merge_n_images=10

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

      data_in = imageset.get_raw_data(i_in)
      assert len(data_in) == 1
      data_in = data_in[0]
      if data_out is None:
        data_out = data_in
      else:
        # FIXME only add pixels to this which are > 0; image pixels < 0
        # are meaningful and should be preserved;
        # Achieved by setting -ve values here to 0 before +=
        # This assumes that -ve values are constant over all images
        data_special = data_in < 0
        data_in.set_selected(data_special, 0)
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
      rsplit = record.split(' ')
      if 'X-Binary-Size:' in record:
        old_size = int(record.split()[-1])
        new_header.append('X-Binary-Size: %d\r\n' % len(compressed))
      elif 'Content-MD5' in record:
        pass
      elif len(rsplit) > 3 and rsplit[1] in { \
        'Exposure_time', 'Angle_increment', 'Exposure_period', 'Count_cutoff', \
        'Phi_increment', 'Omega_increment', 'Chi_increment' }:

        if rsplit[1] == 'Count_cutoff': # needs to be an integer
          new_header.append('%s\n' % ' '.join(rsplit[:2] + ['%d' % (n_images * int(rsplit[2]))] + rsplit[3:]))
        else:
          new_header.append('%s\n' % ' '.join(rsplit[:2] + ['%f' % (n_images * float(rsplit[2]))] + rsplit[3:]))

      else:
        new_header.append('%s\n' % record)

    loop_lines = [n for n, record in enumerate(new_header) if record.startswith('loop_')]
    multiply_fields = { '_diffrn_scan_axis.angle_range', '_diffrn_scan_axis.angle_increment', \
                        '_diffrn_scan_axis.displacement_range', '_diffrn_scan_axis.displacement_increment', \
                        '_diffrn_scan_frame.integration_time', '_diffrn_scan_frame.exposure_time', \
                        '_array_intensities.overload' }
    for loop_start in loop_lines:
      n = loop_start
      modifiers = []
      while True:
        n = n + 1
        line = new_header[n].strip()
        if line in { '', ';' }: # end of loop
          break
        elif line.startswith('_'): # loop header
          if line in multiply_fields:
            modifiers.append(n_images)
          else:
            modifiers.append(None)
        elif any(modifiers): # loop body
          # NOTE: This can break when fields are modified in loops with
          #   'Strings with spaces, as they are seen as multiple columns, or with'
          #   _multiple _columns _defined _on _same _line _they _are _seen _as _one _column
          new_line = [ element if modifier is None else '%f' % (float(element) * modifier) \
                       for modifier, element in zip(modifiers, line.split()) ]
          new_header[n] = '%s\r\n' % ' '.join(new_line)

    tailer = data[data_offset + 4 + old_size:]

    with open(out_image, 'wb') as f:
      f.write(''.join(new_header) + start_tag + compressed + tailer)
      print '%s written' % out_image

  return

def run():
  import libtbx.load_env

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks

  usage = "%s [options] image_*.cbf" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_datablocks_from_images=True,
    epilog=help_message
  )

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)

  n_images = params.merge_n_images
  out_prefix = params.output.image_prefix
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    return

  if len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]

  merge_cbf(imageset, n_images, out_prefix=out_prefix)

if __name__ == '__main__':
  run()
