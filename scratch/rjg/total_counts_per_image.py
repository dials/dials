from __future__ import division

def run(args):
  from dials.array_family import flex
  from dials.util.command_line import Importer
  importer = Importer(args, check_format=True)
  assert len(importer.datablocks) == 1
  imageset = importer.datablocks[0].extract_imagesets()[0]

  total_counts = flex.int()
  for im in imageset:
    total_counts.append(flex.sum(im))

  with open("total_counts_per_image.txt", "wb") as f:
    for i, count in enumerate(total_counts):
      print >> f, "%i %i" %(i, count)

  from matplotlib import pyplot
  pyplot.plot(range(total_counts.size()), total_counts)
  pyplot.xlabel("Image number")
  pyplot.ylabel("Total counts")
  pyplot.show()

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
