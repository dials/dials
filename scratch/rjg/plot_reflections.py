from __future__ import division

def run(args):
  from matplotlib import pyplot
  from dials.util.command_line import Importer
  from scitbx.array_family import flex
  importer = Importer(args)
  sweeps = importer.imagesets
  fig = pyplot.figure()
  pyplot.axes().set_aspect('equal')
  observed_xy = flex.vec2_double()
  predicted_xy = flex.vec2_double()
  for reflection_list in importer.reflections:
    for refl in reflection_list:
      centroid = refl.centroid_position
      if centroid != (0,0,0):
        # this assumes the centroids are given in pixel coordinates
        x, y = centroid[:2]
        observed_xy.append((x,y))
      if refl.image_coord_px != (0, 0):
        x, y = refl.image_coord_px
        predicted_xy.append((x,y))
  obs_x, obs_y = observed_xy.parts()
  pred_x, pred_y = predicted_xy.parts()
  pyplot.scatter(obs_x, obs_y, marker='o', c='red')
  pyplot.scatter(pred_x, pred_y, marker='+', c='blue')
  pyplot.gca().invert_yaxis()
  pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
