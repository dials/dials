from __future__ import division

def calculate_threshold(image, trusted_range):
  from scipy.ndimage.measurements import histogram
  from thresholding import maximum_deviation
  import numpy

  # Cap pixels to within trusted range
  image.shape = -1
  ind = numpy.where(image < trusted_range[0])
  image[ind] = trusted_range[0]
  ind = numpy.where(image > trusted_range[1])
  image[ind] = trusted_range[1]

  # Histogram the pixels
  histo = histogram(image, trusted_range[0], trusted_range[1], trusted_range[1])
  histo = histo / numpy.sum(histo)

  # Calculate the threshold and add to list
  return maximum_deviation(histo)

def select_strong_pixels(sweep, trusted_range):
  from dials.util.command_line import ProgressBar
  import numpy

  # Calculate the threshold
  coordinate = []
  intensity = []
  progress = ProgressBar()
  for i, flex_image in enumerate(sweep):
    image = flex_image.as_numpy_array()
    height, width = image.shape
    threshold = calculate_threshold(image, trusted_range)
    image.shape = -1
    mask = image >= threshold

    ind = numpy.where(mask != 0)[0]
    z = [i] * len(ind)
    y = list(ind / width)
    x = list(ind % width)
    coords = zip(x, y, z)
    coordinate.extend(coords)
    intensity.extend(list(image[ind]))
    progress.update(100.0 * float(i) / len(sweep))
  progress.finished()

  return coordinate, intensity

def find_pixel_nearest_neighbour(pixels, reflections):
  from annlib_ext import AnnAdaptor
  from scitbx.array_family import flex
  from math import sqrt
  import numpy

  # Get the predicted coords
  pred_xyz = []
  for r in reflections:
    x = r.image_coord_px[0]# * 0.172
    y = r.image_coord_px[1]# * 0.172
    z = r.frame_number# * 0.2
    pred_xyz.append((x, y, z))

  # Create the KD Tree
  ann = AnnAdaptor(flex.double(pred_xyz).as_1d(), 3)

  ann.query(flex.double(pixels).as_1d())

#    for i in xrange(len(ann.nn)):
#        print "Neighbor of {0}, index {1} distance {2}".format(
#        obs_xyz[i], ann.nn[i], sqrt(ann.distances[i]))

  return ann.nn.as_numpy_array()

def find_bounding_boxes(index, owner, coordinate):

  width = 9999
  height = 9999
  nframe = 9999
  bbox = [[width, -1, height, -1, nframe, -1]] * (max(owner)+1)

  for r, (x, y, z) in zip(owner, coordinate):
    if x < bbox[r][0]:
      bbox[r][0] = x
    if x + 1 > bbox[r][1]:
      bbox[r][1] = x + 1
    if y < bbox[r][2]:
      bbox[r][2] = y
    if y + 1 > bbox[r][3]:
      bbox[r][3] = y + 1
    if z < bbox[r][4]:
      bbox[r][4] = z
    if z + 1 < bbox[r][5]:
      bbox[r][5] = z + 1

  return bbox


def direction_var(values, weights):
  import numpy
  from scitbx import matrix
  weights = numpy.array(weights)
  valx = numpy.array([x for x, y, z in values])
  valy = numpy.array([y for x, y, z in values])
  valz = numpy.array([z for x, y, z in values])

  # Calculate avergae x, y, z
  avrx = numpy.average(valx, weights=weights)
  avry = numpy.average(valy, weights=weights)
  avrz = numpy.average(valz, weights=weights)

  # Calculate mean direction vector
  s1m = matrix.col((avrx, avry, avrz)).normalize()

  # Calculate angles between vectors
  angles = []
  for s in values:
    angles.append(s1m.angle(s))

  # Calculate variance of angles
  angles = numpy.array(angles)
  var = numpy.dot(weights, (angles)**2)/numpy.sum(weights)
  return var

def find_centroid_and_variance(index, owner, coordinate, intensity, detector):
  from numpy import zeros, int32, argmax, where, average
  from scitbx.array_family import flex
  from scitbx import matrix
  import numpy
  from collections import Counter, defaultdict

  num = Counter(owner)
  groups = defaultdict(list)
  for i, o in enumerate(owner):
    groups[o].append(i)


  var = []
  cent = []
  index = []
  for g in groups:

    xs = []
    ys = []
    zs = []
    s1s = []
    weights = []

    if len(groups) > 6:

      for i in groups[g]:
        x, y, z = coordinate[i]
        c = intensity[i]
        s1 = matrix.col(detector.get_pixel_lab_coord((x+0.5, y+0.5)))
        xs.append(x+0.5)
        ys.append(y+0.5)
        zs.append(z+0.5)
        s1s.append(s1)
        weights.append(c)

      v = direction_var(s1s, weights)
      var.append(v)
      avrx = average(xs, weights=weights)
      avry = average(ys, weights=weights)
      avrz = average(zs, weights=weights)
      cent.append((avrx, avry, avrz))
      index.append(g)

  # Return a list of centroids and variances
  return index, numpy.array(cent), numpy.array(var)


def filter_objects_by_distance(index, centroid, reflections, max_distance):

  import numpy
  from math import sqrt

  # Only accept objects closer than max distance
  indices = []
  for i, (r, (xc, yc, zc)) in enumerate(zip(index, centroid)):
    x, y = reflections[int(r)].image_coord_px
    z = reflections[int(r)].frame_number
    if sqrt((x - xc)**2 + (y - yc)**2 + (z - zc)**2) > max_distance:
      continue
    else:
      indices.append(i)

  return numpy.array(indices)


def calculate_sigma_beam_divergence(var):
  '''Calculate the beam divergence as the sum of centroid variance of the
  intensity weighted diffracted beam directions.'''
  from math import sqrt

  # Calculate the sum of s^2
  sum_variance = reduce(lambda x, y: x + y, var)

  # Return the beam divergence as the sum / num reflections
  return sqrt(sum_variance / len(var))


if __name__ == '__main__':

  import os
  import libtbx.load_env
  from dials.util.nexus import NexusFile
  from glob import glob
  from dxtbx.sweep import SweepFactory
  from math import pi
  from scitbx import matrix

  try:
    dials_regression = libtbx.env.dist_path('dials_regression')
  except KeyError, e:
    print 'FAIL: dials_regression not configured'
    raise

  # The XDS values
  xds_sigma_d = 0.060
  xds_sigma_m = 0.154

  # Read the reflections file
  filename = os.path.join(dials_regression,
      'centroid_test_data', 'reflections.h5')
  #filename = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'reflections.h5')
  handle = NexusFile(filename, 'r')

  # Get the reflection list
  print 'Reading reflections.'
  reflections = handle.get_reflections()
  print 'Read {0} reflections.'.format(len(reflections))

  # Read images
  template = os.path.join(dials_regression,
      'centroid_test_data', 'centroid_*.cbf')
  #template = os.path.join('/home/upc86896/Data/X1_strong_M1S1_1_', 'X1_strong_M1S1_1_*.cbf')
  filenames = glob(template)

  # Load the sweep
  print 'Loading sweep'
  sweep = SweepFactory.sweep(filenames)
  print 'Loaded sweep of {0} images.'.format(len(sweep))

  # Select the strong pixels to use in the divergence calculation
  print 'Select the strong pixels from the images.'
  trusted_range = (0, 20000)
  coordinate, intensity = select_strong_pixels(sweep, trusted_range)
  print 'Selected {0} pixels'.format(len(coordinate))

  print 'Finding pixel nearest neighbours.'
  owner = find_pixel_nearest_neighbour(coordinate, reflections)
  print 'Found {0} nearest neighbours.'.format(len(owner))

  print 'Grouping pixels by reflection.'
  index = range(len(owner))
  index = sorted(index, key=lambda x: owner[x])
  print 'Sorted {0} pixels.'.format(len(index))

  print 'Finding bounding boxes.'
  bbox = find_bounding_boxes(index, owner, coordinate)
  print 'Found {0} bounding boxes'.format(len(bbox))

  print 'Find centroid and variance'
  index, centroid, variance = find_centroid_and_variance(index, owner, coordinate,
      intensity, sweep.get_detector())

  print 'Filter objects by distance from nearest reflection.'
  max_distance = 2
  indices = filter_objects_by_distance(index, centroid, reflections, max_distance)
  print '{0} remaining objects'.format(len(indices))

  print 'Calculate the e.s.d of the beam divergence'
  sigma_d = calculate_sigma_beam_divergence(variance[indices])
  print 'Sigma_d = {0} deg'.format(sigma_d * 180.0 / pi)
  print 'XDS Sigma_d = {0} deg'.format(xds_sigma_d)
