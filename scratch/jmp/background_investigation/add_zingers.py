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

  dials.python add_zingers.py image_*.cbf

'''

phil_scope = iotbx.phil.parse("""\

zingers {

  average_per_image = 10
    .type = float(value_min=0)
    .help = "Average zingers per image (poisson mean)"

  average_intensity = 1000
    .type = float(value_min=0)
    .help = "Average intensity of zingers (poisson mean)"

}

output {
  directory = "."
    .type = str
    .help = "The output directory"
}
""", process_includes=True)


def add_zingers(imageset, params):

  from dxtbx.format.FormatCBF import FormatCBF
  assert issubclass(imageset.reader().get_format_class(), FormatCBF), (
    "Only CBF format images supported")

  from cbflib_adaptbx import compress
  import binascii
  from numpy.random import poisson
  from random import sample
  import os.path

  for i in range(len(imageset)):
    image_data = imageset[i]

    num = poisson(params.zingers.average_per_image)
    index = sample(range(len(image_data)), num)
    value = list(poisson(params.zingers.average_intensity, num))
    for j, v in zip(index, value):
      image_data[j] += v
    out_image = os.path.join(params.output.directory, "image_%04i.cbf" % i)

    start_tag = binascii.unhexlify('0c1a04d5')

    data = open(imageset.get_path(i), 'rb').read()
    data_offset = data.find(start_tag)
    cbf_header = data[:data_offset]

    new_header = []
    compressed = compress(image_data)

    old_size = 0

    for record in cbf_header.split('\n')[:-1]:
      if 'X-Binary-Size:' in record:
        old_size = int(record.split()[-1])
        new_header.append('X-Binary-Size: %d\r\n' % len(compressed))
      elif 'Content-MD5' in record:
        pass
      else:
        new_header.append('%s\n' % record)

    tailer = data[data_offset + 4 + old_size:]

    with open(out_image, 'wb') as f:
      f.write(''.join(new_header) + start_tag + compressed + tailer)
      print '%s written with %d zingers' % (out_image, num)

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

  params, options = parser.parse_args(show_diff_phil=True)
  datablocks = flatten_datablocks(params.input.datablock)
  if len(datablocks) > 1:
    raise Sorry("Only one DataBlock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]

  try:
    import os
    os.makedirs(params.output.directory)
  except Exception:
    pass

  add_zingers(imageset, params)
