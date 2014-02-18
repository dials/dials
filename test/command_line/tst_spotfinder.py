from __future__ import division
import os
import cPickle as pickle
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from dials.model.data import ReflectionList # import dependency

def exercise_spotfinder():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_spotfinder: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "centroid*.cbf")
  args = ["dials.find_spots", template, "-o spotfinder.pickle", "--nproc=1",
          "save_shoeboxes=True"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 459
    refl = reflections[0]
    assert approx_equal(refl['intensity.raw.value'], 142)
    assert approx_equal(refl['bbox'], (1258, 1260, 537, 541, 0, 1))
    assert approx_equal(refl['xyzobs.px.value'],
                        (1258.7957746478874, 539.112676056338, 0.5))
    assert "shoebox" in reflections
  print 'OK'

  # now with a resolution filter
  args = ["dials.find_spots", "d_min=2", "d_max=15",
          template, "-o spotfinder.pickle", "--nproc=1", "save_shoeboxes=False"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 371
    refl = reflections[0]
    assert "shoebox" not in reflections
  print 'OK'

  # now with more generous parameters
  args = ["dials.find_spots", "min_spot_size=3", "max_separation=3",
          template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 679
  print 'OK'

  # now with XFEL stills
  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/spotfinding_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "idx-s00-20131106040302615.cbf")
  args = ["dials.find_spots", template, "-o spotfinder.pickle", "--nproc=1"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 2643
  print 'OK'

def exercise_polygon():
  from dials.algorithms.peak_finding.spotfinder_factory import polygon

  x = 1
  y = 1
  vertices = [(0,0), (2,0), (2,2), (0,2)]
  poly = polygon(vertices)
  assert poly.is_inside(x,y)

  poly = polygon([(3,5), (40,90), (80,70), (50,50), (70,20)])
  for p in [(42,40), (16,25), (64,67), (16,30), (48, 45), (21,30)]:
    assert poly.is_inside(p[0], p[1])
  for p in [(59,4), (70,15), (57,14), (21,78), (37,100), (88,89)]:
    assert not poly.is_inside(p[0], p[1])

  if 0:
    # for visual confirmation of algorithm
    from scitbx.array_family import flex
    inside_points = flex.vec2_double()
    outside_points = flex.vec2_double()
    import random
    x_max = 100
    y_max = 100
    for i in range(1000):
      x = random.randint(0,x_max)
      y = random.randint(0,y_max)
      is_inside = poly.is_inside(x, y)
      if is_inside:
        inside_points.append((x,y))
      else:
        outside_points.append((x,y))
    from matplotlib import pyplot
    from matplotlib.patches import Polygon
    v = poly.vertices + poly.vertices[:1]
    fig = pyplot.figure()
    ax = fig.add_subplot(111)
    ax.add_patch(Polygon(poly.vertices, closed=True, fill=False))
    inside_x, inside_y = inside_points.parts()
    outside_x, outside_y = outside_points.parts()
    ax.scatter(inside_x, inside_y, marker='+', c='r')
    ax.scatter(outside_x, outside_y, marker='+', c='b')
    pyplot.show()


def run():
  exercise_polygon()
  exercise_spotfinder()

if __name__ == '__main__':
  run()
  print 'OK'
