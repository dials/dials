#!/usr/bin/env python
# convert_to_cbf.py
#
#   Copyright (C) 2018 Diamond Light Source, Graeme Winter
#
#   This code is distributed under the BSD license, a copy of which is
#   included in the root directory of this package.
#

from __future__ import absolute_import, division, print_function

import iotbx.phil

help_message = '''

Convert data which can be read by DIALS, given a datablock, to CBF format -
with ersatz miniCBF header. Can be used e.g. with HDF5 format data.

Examples::

  dials.convert_to_cbf datablock.json prefix=data_as_cbf

'''

phil_scope = iotbx.phil.parse("""\
output {
  template = as_cbf_%04d.cbf
    .type = path
}
""", process_includes=True)

def convert_to_cbf(imageset, template):
  from cbflib_adaptbx import compress
  import binascii
  start_tag = binascii.unhexlify('0c1a04d5')

  for i in range(len(imageset)):
    data = imageset.get_raw_data(i)[0]
    compressed = compress(data)

    header = '''###CBF: VERSION 3.14 from DIALS conversion

data_data_%06d

_array_data.header_convention "SLS_1.0"
_array_data.header_contents
;
# WARNING: FOR XDS PROCESSING ONLY. MAY EAT YOUR KITTEN. YOU HAVE BEEN WARNED
;

_array_data.data
;
--CIF-BINARY-FORMAT-SECTION--
Content-Type: application/octet-stream;
     conversions="x-CBF_BYTE_OFFSET"
Content-Transfer-Encoding: BINARY
X-Binary-Size: %d
X-Binary-ID: 1
X-Binary-Element-Type: "signed 32-bit integer"
X-Binary-Element-Byte-Order: LITTLE_ENDIAN
X-Binary-Number-of-Elements: %d
X-Binary-Size-Fastest-Dimension: %d
X-Binary-Size-Second-Dimension: %d
X-Binary-Size-Padding: 0

''' % (i, len(compressed), data.size(), data.focus()[1], data.focus()[0])

    with open(template % (i + 1), 'wb') as f:
      print(template % (i + 1))
      f.write(''.join(header) + start_tag + compressed)

  return

def run():
  import libtbx.load_env

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks

  usage = "%s [options] datablock.json" % libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    epilog=help_message
  )

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)

  template = params.output.template
  datablocks = flatten_datablocks(params.input.datablock)

  if len(datablocks) == 0:
    parser.print_help()
    return

  if len(datablocks) > 1:
    raise Sorry("Only one datablock can be processed at a time")
  else:
    imagesets = datablocks[0].extract_imagesets()
    assert len(imagesets) == 1
    imageset = imagesets[0]

  convert_to_cbf(imageset, template)

if __name__ == '__main__':
  run()
