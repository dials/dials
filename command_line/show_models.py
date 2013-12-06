from __future__ import division

def run(args):
  from dials.util.command_line import Importer
  importer = Importer(args)
  for crystal_model in importer.crystals:
    print crystal_model
  for imageset in importer.imagesets:
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
