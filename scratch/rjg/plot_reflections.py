from __future__ import division

import iotbx.phil

master_phil_scope = iotbx.phil.parse("""
scan_range = None
  .help = "The range of images to use in indexing. Number of arguments"
          "must be a factor of two. Specifying \"0 0\" will use all images"
          "by default. The given range follows C conventions"
          "(e.g. j0 <= j < j1)."
  .type = ints(size=2)
  .multiple = True
""")

def run(args):
  from matplotlib import pyplot
  from dials.util.command_line import Importer
  from scitbx.array_family import flex
  from libtbx.phil import command_line
  importer = Importer(args)
  assert len(importer.imagesets) == 1
  imageset = importer.imagesets[0]
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  args = importer.unhandled_arguments
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  params = working_phil.extract()
  fig = pyplot.figure()
  pyplot.axes().set_aspect('equal')
  observed_xy = flex.vec2_double()
  predicted_xy = flex.vec2_double()
  for reflection_list in importer.reflections:
    for refl in reflection_list:
      if len(params.scan_range):
        reflections_in_range = False
        for scan_range in params.scan_range:
          if scan_range is None: continue
          range_start, range_end = scan_range
          if refl.frame_number >= range_start and refl.frame_number < range_end:
            reflections_in_range = True
            break
        if not reflections_in_range:
          continue
      centroid_position = refl.centroid_position
      centroid_variance = refl.centroid_variance
      if centroid_position != (0,0,0):
        # just a quick check for now that the reflections haven't come from
        # somewhere else
        assert refl.image_coord_mm == (0,0)
        # this assumes the centroids are given in pixel coordinates
        from dials.algorithms.centroid import centroid_px_to_mm_panel
        centroid_position, centroid_variance, _ = centroid_px_to_mm_panel(
          detector[refl.panel_number], scan,
          centroid_position,
          centroid_variance,
          (1,1,1))
        x, y = centroid_position[:2]
        observed_xy.append((x,y))
      if refl.image_coord_px != (0, 0):
        x, y = refl.image_coord_px
        predicted_xy.append((x,y))
  obs_x, obs_y = observed_xy.parts()
  pred_x, pred_y = predicted_xy.parts()
  pyplot.scatter(obs_x, obs_y, marker='o', c='red')
  pyplot.scatter(pred_x, pred_y, marker='+', c='blue')
  assert len(detector) == 1
  panel = detector[0]
  pyplot.xlim(0, panel.get_image_size_mm()[0])
  pyplot.ylim(0, panel.get_image_size_mm()[1])
  pyplot.gca().invert_yaxis()
  pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
