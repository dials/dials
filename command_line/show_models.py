from __future__ import division

def run(args):


  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks

  parser = OptionParser(
    read_experiments=True,
    read_datablocks=True,
    check_format=False)

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
      crystal.show(show_scan_varying=True)
  if datablocks is not None:
    for datablock in datablocks:
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
  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
