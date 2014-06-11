from __future__ import division

def run(args):
  from dials.util.command_line import Importer
  importer = Importer(args)
  experiments = importer.experiments
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
      print crystal
  if importer.datablocks is not None:
    for datablock in importer.datablocks:
      imagesets = datablock.extract_imagesets()
      for imageset in imagesets:
        print imageset
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
