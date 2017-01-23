from __future__ import absolute_import, division
import os
from dials.array_family import flex # import dependency
import cPickle as pickle
import libtbx.load_env
from libtbx import easy_run
from libtbx.test_utils import approx_equal
from glob import glob

def exercise_spotfinder():
  if not libtbx.env.has_module("dials_regression"):
    print "Skipping exercise_spotfinder: dials_regression not present"
    return

  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/centroid_test_data",
    test=os.path.isdir)
  template = glob(os.path.join(data_dir, "centroid*.cbf"))
  args = ["dials.find_spots", ' '.join(template), "output.reflections=spotfinder.pickle",
          "output.shoeboxes=True"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 653, len(reflections)
    refl = reflections[0]
    assert approx_equal(refl['intensity.sum.value'], 42)
    assert approx_equal(refl['bbox'], (1398, 1400, 513, 515, 0, 1))
    assert approx_equal(refl['xyzobs.px.value'],
                        (1399.1190476190477, 514.2142857142857, 0.5))
    assert "shoebox" in reflections
  print 'OK'

  # now with a resolution filter
  args = ["dials.find_spots", "filter.d_min=2", "filter.d_max=15",
          ' '.join(template), "output.reflections=spotfinder.pickle", "output.shoeboxes=False"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 467, len(reflections)
    assert "shoebox" not in reflections
  print 'OK'

  # now with more generous parameters
  args = ["dials.find_spots", "min_spot_size=3",
          "max_separation=3",
          ' '.join(template), "output.reflections=spotfinder.pickle"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 678, len(reflections)
  print 'OK'

  # Now with a user defined mask
  template = glob(os.path.join(data_dir, "centroid*.cbf"))
  args = ["dials.find_spots", ' '.join(template), "output.reflections=spotfinder.pickle",
          "output.shoeboxes=True",
          "lookup.mask=%s" % os.path.join(data_dir, "mask.pickle")]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_json_file(os.path.join(data_dir,
                                                              "datablock.json"))
    assert(len(datablocks) == 1)
    imageset = datablocks[0].extract_imagesets()[0]
    detector = imageset.get_detector()
    beam = imageset.get_beam()
    for x, y, z in reflections['xyzobs.px.value']:
      d = detector[0].get_resolution_at_pixel(beam.get_s0(), (x, y))
      assert(d >= 3)

  # Now with a user defined mask
  template = glob(os.path.join(data_dir, "centroid*.cbf"))
  args = ["dials.find_spots", ' '.join(template), "output.reflections=spotfinder.pickle",
          "output.shoeboxes=True",
          "region_of_interest=800,1200,800,1200"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    x, y, z = reflections['xyzobs.px.value'].parts()
    assert x.all_ge(800)
    assert y.all_ge(800)
    assert x.all_lt(1200)
    assert y.all_lt(1200)

  print 'OK'


  # now with XFEL stills
  data_dir = libtbx.env.find_in_repositories(
    relative_path="dials_regression/spotfinding_test_data",
    test=os.path.isdir)
  template = os.path.join(data_dir, "idx-s00-20131106040302615.cbf")
  args = ["dials.find_spots", template, "output.reflections=spotfinder.pickle"]
  result = easy_run.fully_buffered(command=" ".join(args)).raise_if_errors()
  assert os.path.exists("spotfinder.pickle")
  with open("spotfinder.pickle", "rb") as f:
    reflections = pickle.load(f)
    assert len(reflections) == 2643, len(reflections)
  print 'OK'



def exercise_polygon():
  from dials.algorithms.polygon import polygon

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
  from dials.test import cd_auto
  with cd_auto(__file__):
    run()
    print 'OK'
