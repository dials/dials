from __future__ import division

def run(args):
  from dials.util.command_line import Importer
  importer = Importer(args)
  for crystal_model in importer.crystals:
    print crystal_model
  for imageset in importer.imagesets:
    print imageset
  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
