from __future__ import division

if __name__ == '__main__':

  from dials.util.command_line import Importer
  from dials.array_family import flex
  import sys

  # Get the imageset
  importer = Importer(sys.argv[1:])
  assert(len(importer.datablocks) == 1)
  imagesets = importer.datablocks[0].extract_imagesets()
  assert(len(imagesets) == 1)
  imageset = imagesets[0]

  # Create a mask
  image = imageset[0]
  mask = image >= 0

  # Loop through all the images and sum the pixels
  counts = []
  for i, image in enumerate(imageset):
    print 'Processing image %d' % i
    counts.append(flex.sum(image.select(mask.as_1d())))

  # Write the counts to file
  with open("counts_per_image.txt", "w") as outfile:
    for item in counts:
      outfile.write('%s\n' % str(item))

  from matplotlib import pylab
  pylab.plot(counts)
  pylab.show()
