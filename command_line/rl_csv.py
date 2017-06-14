# LIBTBX_SET_DISPATCHER_NAME dev.dials.csv
from __future__ import absolute_import, division

import iotbx.phil
from scitbx.array_family import flex
from scitbx import matrix
from dials.util.options import OptionParser
from dials.util.options import flatten_datablocks, flatten_reflections
from dials.algorithms.indexing.indexer \
     import indexer_base, filter_reflections_by_scan_range

import libtbx.load_env

phil_scope = iotbx.phil.parse("""
output {
  csv = rl.csv
    .type = path
    .help = 'Output filename for reciprocal mapped reflections'
}
""")

master_params = phil_scope.fetch().extract()

help_message = "%s datablock.json strong.pickle output.csv=rl.csv" % \
  libtbx.env.dispatcher_name

def run(args):
  import libtbx.load_env
  from dials.util import log
  usage = "%s [options] datablock.json strong.pickle" %libtbx.env.dispatcher_name

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_datablocks=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0 or len(reflections) == 0:
    parser.print_help()
    exit(0)

  imagesets = []

  for db in datablocks:
    imagesets.extend(db.extract_imagesets())

  assert len(imagesets) == len(reflections)

  fout = open(params.output.csv, 'w')
  fout.write('# x,y,z,experiment_id,imageset_id\n')

  for k, (imageset, refl) in enumerate(zip(imagesets, reflections)):
    if 'imageset_id' not in refl:
      refl['imageset_id'] = refl['id']

    reflmm = indexer_base.map_spots_pixel_to_mm_rad(
      spots=refl, detector=imageset.get_detector(), scan=imageset.get_scan())

    indexer_base.map_centroids_to_reciprocal_space(
      reflmm, detector=imageset.get_detector(), beam=imageset.get_beam(),
      goniometer=imageset.get_goniometer())

    rlp = reflmm['rlp']

    for _rlp in rlp:
      fout.write('%f,%f,%f,%d,%d\n' % (_rlp[0], _rlp[1], _rlp[2], k, k))

    print 'Appended %d reflections to %s' % (len(rlp), params.output.csv)

  fout.close()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
