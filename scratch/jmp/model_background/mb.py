
from __future__ import division


# def online_variance(data):
#     n = 0
#     mean = 0.0
#     M2 = 0.0

#     for x in data:
#         n = n + 1
#         delta = x - mean
#         mean = mean + delta/n
#         M2 = M2 + delta*(x - mean)

#     if n < 2:
#         return float('nan');
#     else:
#         return M2 / (n - 1)

# class RunningMeanAndVariance(object):

#   def __init__(self, size):
#     from dials.array_family import flex
#     self.n = flex.int(flex.grid(size), 0)
#     self.mean = flex.double(flex.grid(size), 0.0)
#     self.m2 = flex.double(flex.grid(size), 0.0)

#   def add(self, image, mask):
#     m = (mask == True).as_1d().as_int()
#     x = data.as_double() * m.as_double()
#     self.n += m
#     delta = x - self.mean
#     n1 = flex.double(image.accessor(), 0)
#     self.mean = self.mean + delta * n1
#     self.m2 = self.m2 + delta * (x - mean)

#   def variance(self):
#     return self.m2 / (self.n - 1)

#   def stddev(self):
#     from dials.array_family import flex
#     return flex.sqrt(self.variance())


if __name__ == '__main__':

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.array_family import flex
  from collections import defaultdict
  from dials.algorithms.shoebox import MaskCode

  # Create the option parser
  parser = OptionParser(
    read_experiments=True,
    read_reflections=True)

  # Parse the arguments
  params, options = parser.parse_args(show_diff_phil=True)

  # Get the experiments and reflections
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  # Get the imageset
  assert len(experiments) == 1
  assert len(reflections) == 1
  reflections = reflections[0]
  imageset = experiments[0].imageset

  # Create the reflection lookup
  bbox = reflections['bbox']
  reflection_lookup = defaultdict(list)
  for i in range(len(bbox)):
    for j in range(bbox[i][4],bbox[i][5]):
      reflection_lookup[j].append(i)

  width, height = experiments[0].detector[0].get_image_size()
  sum_background = flex.double(flex.grid(height, width), 0)
  sum_sq_background = flex.double(flex.grid(height, width), 0)
  count = flex.int(flex.grid(height, width), 0)

  # Loop through all images
  print "START"
  for frame in range(len(imageset)):

    # Get the subset of reflections on this image and compute the mask
    subset = reflections.select(flex.size_t(reflection_lookup[frame]))

    subset['shoebox'] = flex.shoebox(
      subset['panel'],
      subset['bbox'],
      allocate=True)

    subset.compute_mask(experiments)

    # Get the mask and data
    mask = imageset.get_mask(frame)[0]
    data = imageset.get_raw_data(frame)[0]

    sbox_mask = subset['shoebox'].apply_background_mask(frame, 1, (height, width))

    mask = mask & sbox_mask

    #from dials.algorithms.image.threshold import DispersionThreshold
    #threshold = DispersionThreshold(
    #  data.all(),
    #  (3,3),
    #  6,3,0,2)
    #new_mask = flex.bool(mask.accessor())
    #threshold(data, mask, new_mask)
    #mask = mask & (~new_mask)

    import cPickle as pickle
    pickle.dump((data, mask), open("first_image.pickle", "w"))
    exit(0)
    m = (mask == True).as_1d().as_int()
    x = data.as_double() * m.as_double()
    sum_background += x
    sum_sq_background += x * x
    count += m

    average = flex.sum(sum_background) / flex.sum(count)

    print "Image %d: selected %d reflections, avr=%f" % (
      frame,
      len(subset),
      average)

    # from matplotlib import pylab
    # pylab.imshow((count > 0).as_numpy_array())
    # pylab.show()

  average = flex.double(len(sum_background))
  variance = flex.double(len(sum_background))
  count_mask = count > 1
  indices = flex.size_t(range(len(mask))).select(count_mask.as_1d())
  from matplotlib import pylab
  pylab.imshow(count_mask.as_numpy_array())
  pylab.show()

  sumb = sum_background.as_1d().select(indices)
  numb = count.as_1d().select(indices).as_double()
  avrb = sumb / numb
  sumsqb = sum_sq_background.as_1d().select(indices)
  varb = (sumsqb - sumb*sumb / numb) / (numb - 1)
  average.set_selected(indices, avrb)
  average.reshape(count_mask.accessor())
  variance.set_selected(indices, varb)
  variance.reshape(count_mask.accessor())

  print "Saving to model.pickle"
  with open("model.pickle", "w") as outfile:
    import cPickle as pickle
    pickle.dump((average, count_mask, flex.sqrt(variance)), outfile)
