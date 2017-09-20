from __future__ import absolute_import, division

def evaluate_gaussian(x, a, x0, sx):

  from math import exp

  assert(len(x) == len(x0))
  assert(len(x) == len(sx))

  g = 0.0
  for xi, x0i, sxi in zip(x, x0, sx):
    g += (xi - x0i)**2 / (2.0 * sxi**2)

  return a * exp(-g)

def gaussian(size, a, x0, sx):

  from scitbx.array_family import flex

  result = flex.double(flex.grid(size))

  index = [0] * len(size)
  while True:
    result[index] = evaluate_gaussian(index, a, x0, sx)
    for j in range(len(size)):
      index[j] += 1
      if index[j] < size[j]:
        break
      index[j] = 0
      if j == len(size) - 1:
        return result
