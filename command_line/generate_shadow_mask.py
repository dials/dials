# LIBTBX_SET_DISPATCHER_NAME dev.dials.generate_shadow_mask
from __future__ import absolute_import, division, print_function
from scitbx.array_family import flex

def main():
  import sys
  run(sys.argv[1:])

def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_datablocks, flatten_experiments
  import libtbx.load_env

  usage = "%s [options] (datablock.json|experiments.json)" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    read_datablocks=True,
    read_datablocks_from_images=True,
    read_experiments=True)

  params, options = parser.parse_args()

  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)
  assert[len(datablocks), len(experiments)].count(1) == 1

  if datablocks:
    datablock = datablocks[0]
    imagesets = datablock.extract_imagesets()
  else:
    imagesets = experiments.imagesets()

  imageset = imagesets[0]
  goniometer = imageset.get_goniometer()
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  format_class = imageset.masker().format_class(imageset.paths()[0])
  masker = format_class.get_goniometer_shadow_masker()
  angles = goniometer.get_angles()
  start, end = scan.get_oscillation_range()
  step = scan.get_oscillation()[1]

  scan_points = flex.double(libtbx.utils.frange(start, end, step=step))
  array_range = imageset.get_array_range()
  image_size = detector[0].get_image_size()

  nz = array_range[1] - array_range[0]
  nx, ny = image_size

  import h5py
  import numpy
  fout = h5py.File('shadow.hdf5', 'w')
  shadow = fout.create_dataset('shadow', (nz, ny, nx), fillvalue=0,
                               dtype='i1', chunks=(1, ny, nx),
                               compression='gzip')

  from dials.util.command_line import ProgressBar
  pb = ProgressBar(title="Casting shadows")
  for j, scan_angle in enumerate(scan_points):
    pb.update(100 * j / nz)
    mask = masker.get_mask(detector, scan_angle=scan_angle)
    for p_id in range(len(detector)):
      if mask[p_id]:
        slab = (~mask[p_id]).as_numpy_array().astype(numpy.int8)
        shadow[j,:,:] = slab
  pb.finished()
  fout.close()

if __name__ == '__main__':
  main()
