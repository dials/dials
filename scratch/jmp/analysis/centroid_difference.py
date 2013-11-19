from __future__ import division

class Script(object):

  def __init__(self):
    pass


  def run(self, filename):
    from dials.model.data import ReflectionList, Reflection
    from math import sqrt
    import cPickle as pickle
    import numpy

    rlist = pickle.load(open(filename, "rb"))
    print "read file"

    r_to_use = [r for r in rlist
                    if r.is_valid() and
                       r.bounding_box[1] - r.bounding_box[0] >= 5 and
                       r.bounding_box[3] - r.bounding_box[2] >= 5 and
                       r.intensity / sqrt(r.intensity_variance) > 10]
    print "Filtered list", len(r_to_use)


    px_prd = [r.image_coord_px for r in r_to_use]
    f_prd = [r.frame_number for r in r_to_use]
    px_obs = [r.centroid_position for r in r_to_use]

    x_prd, y_prd = zip(*px_prd)
    x_obs, y_obs, f_obs = zip(*px_obs)

    diff = []
    for xp, yp, zp, xo, yo, zo in zip(x_prd, y_prd, f_prd, x_obs, y_obs, f_obs):
      diff.append(sqrt((xp - xo)**2 + (yp - yo)**2 + (zp - zo)**2))

    x_diff = [p - o for p, o in zip(x_prd, x_obs)]
    y_diff = [p - o for p, o in zip(y_prd, y_obs)]
    f_diff = [p - o for p, o in zip(f_prd, f_obs)]

    from matplotlib import pylab
    pylab.subplot(3, 1, 1)
    pylab.scatter(x_prd, x_diff)
    pylab.axhline(numpy.mean(x_diff))
    pylab.subplot(3, 1, 2)
    pylab.scatter(y_prd, y_diff)
    pylab.axhline(numpy.mean(y_diff))
    pylab.subplot(3, 1, 3)
    pylab.scatter(f_prd, f_diff)
    pylab.axhline(numpy.mean(f_diff))
    pylab.show()

    #from scipy.interpolate import griddata
    points = zip(y_prd, x_prd)
    pickle.dump((points, x_diff, y_diff), open("temp_data.pickle", "wb"))
    #print "gridiing"
    #x_grid = griddata(points, x_diff)
    #print "gridded"
    #y_grid = griddata(points, y_diff)

    #pylab.subplot(2, 1, 1)
    #pylab.imshow(x_grid)
    #pylab.subplot(2, 1, 2)
    #pylab.imshow(y_grid)
    #pylab.show()

def scipy_stuff():
  from scipy.interpolate import griddata
  from matplotlib import pylab
  import cPickle as pickle
  print "loading points"
  points, x_diff, y_diff = pickle.load(open("temp_data.pickle", "rb"))

  y_pts, x_pts = zip(*points)

  print "Creating grid points"
  grid_points = []
  for j in range(2500):
    for i in range(2500):
      grid_points.append((j, i))

  print "Gridding data"
  x_grid = griddata(points, x_diff, grid_points)
  y_grid = griddata(points, y_diff, grid_points)
  x_grid.shape = (2500, 2500)
  y_grid.shape = (2500, 2500)

  print "Plotting"
  pylab.subplot(3, 1, 1)
  pylab.imshow(x_grid)
  pylab.subplot(3, 1, 2)
  pylab.imshow(y_grid)
  pylab.subplot(3, 1, 3)
  pylab.scatter(x_pts, y_pts)
  pylab.show()


if __name__ == '__main__':
  import sys
  script = Script()
  script.run(sys.argv[1])
#    scipy_stuff()
