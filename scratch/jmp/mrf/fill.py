
from __future__ import division

def generate_image(function, height, width):
  from dials.array_family import flex
  image = flex.double(flex.grid(height, width))
  for j in range(height):
    for i in range(width):
      image[j,i] = function(j,i)
  return image

def generate_mask(height, width, y0, y1, x0, x1):
  from dials.array_family import flex
  mask = flex.bool(flex.grid(height, width), True)
  mask[y0:y1,x0:x1] = flex.bool(flex.grid(y1-y0,x1-x0), False)
  return mask

def apply_mask(image, mask):
  mask_double = mask.as_1d().as_double()
  mask_double.reshape(mask.accessor())
  return image * mask_double

def plane(mx, my, c):
  def function(x, y):
    return mx*x + my*y + c
  return function


if __name__ == '__main__':

  from matplotlib import pylab

  height = 20
  width = 20

  # Generate image
  image = generate_image(plane(1, 1, 10), height, width)

  # Set mask
  mask = generate_mask(height, width, 7, 13, 7, 13)
  image = apply_mask(image, mask)

  # Fill pixels

  pylab.imshow(image.as_numpy_array(), interpolation='none')
  pylab.show()
