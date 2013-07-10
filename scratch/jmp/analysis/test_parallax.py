
import sys
from dxtbx.serialize import load
from dxtbx.model import ParallaxCorrectedPxMmStrategy, SimplePxMmStrategy
from scitbx.array_family import flex
from matplotlib import pylab
convert = ParallaxCorrectedPxMmStrategy(0.252500934883)
convert2 = SimplePxMmStrategy()

sweep = load.imageset(sys.argv[1])

detector = sweep.get_detector()

#print detector[0].pixel_to_millimeter((0, 1))

image_size = sweep[0].all()
image = flex.double(flex.grid(image_size))
for j in range(image_size[0]):
    for i in range(image_size[1]):
        mm1 = convert2.to_millimeter(detector[0], (i, j))
        px2 = convert.to_pixel(detector[0], mm1)
        image[j,i] = j - px2[1]


print flex.max(image), flex.min(image)

pylab.imshow(image.as_numpy_array())
pylab.show()
