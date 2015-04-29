from __future__ import division

import iotbx.phil

help_message = '''

Examples::

  dials.show_models datablock.json

  dials.show_models experiments.json

  dials.show_models image_*.cbf

'''

phil_scope = iotbx.phil.parse("""\
show_scan_varying = False
  .type = bool
  .help = "Whether or not to show the crystal at each scan point."
""", process_includes=True)


def run(args):

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
  from libtbx.utils import Abort
  import libtbx.load_env

  usage = "%s [options] datablock.json | experiments.json | image_*.cbf" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_datablocks=True,
    read_datablocks_from_images=True,
    check_format=False,
    epilog=help)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)

  if experiments is not None:
    for detector in experiments.detectors():
      print detector
    for beam in experiments.beams():
      print beam
    for scan in experiments.scans():
      print scan
    for goniometer in experiments.goniometers():
      print goniometer
    for crystal in experiments.crystals():
      crystal.show(show_scan_varying=params.show_scan_varying)
      if crystal.num_scan_points:
        from scitbx.array_family import flex
        from cctbx import uctbx
        abc = flex.vec3_double()
        angles = flex.vec3_double()
        for n in range(crystal.num_scan_points):
          a, b, c, alpha, beta, gamma = crystal.get_unit_cell_at_scan_point(n).parameters()
          abc.append((a, b, c))
          angles.append((alpha, beta, gamma))
        a, b, c = abc.mean()
        alpha, beta, gamma = angles.mean()
        mean_unit_cell = uctbx.unit_cell((a, b, c, alpha, beta, gamma))
        print "  Average unit cell: %s" %mean_unit_cell
  if datablocks is not None:
    for datablock in datablocks:
      if datablock.format_class() is not None:
        print 'Format: %s' %datablock.format_class()
      imagesets = datablock.extract_imagesets()
      for imageset in imagesets:
        try: print imageset.get_template()
        except Exception: pass
        print imageset.get_detector()
        print imageset.get_beam()
        if imageset.get_scan() is not None:
          print imageset.get_scan()
        if imageset.get_goniometer() is not None:
          print imageset.get_goniometer()
  if experiments is None and datablocks is None:
    raise Abort('No experiments or datablocks specified')
  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
