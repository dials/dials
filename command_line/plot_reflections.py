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
first_n_reflections = None
  .type = int(value_min=0)
output {
  file_name = centroids.png
    .type = path
  show_plot = False
    .type = bool
  dpi = 300
    .type = int(value_min=0)
  size_inches = 10,10
    .type = floats(size=2, value_min=0)
}
""")

def run(args):
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
  observed_xy = flex.vec2_double()
  predicted_xy = flex.vec2_double()
  for reflection_list in importer.reflections:
    if len(params.scan_range):
      sel = flex.bool(reflection_list.size(), False)
      centroid_positions = reflection_list.centroid_position()
      centroids_frame = centroid_positions.parts()[2]
      reflections_in_range = False
      for scan_range in params.scan_range:
        if scan_range is None: continue
        range_start, range_end = scan_range
        sel |= ((centroids_frame >= range_start) & (centroids_frame < range_end))
      reflection_list = reflection_list.select(sel)
    if params.first_n_reflections is not None:
      centroid_positions = reflection_list.centroid_position()
      centroids_frame = centroid_positions.parts()[2]
      perm = flex.sort_permutation(centroids_frame)
      perm = perm[:min(reflection_list.size(), params.first_n_reflections)]
      print flex.max(centroids_frame.select(perm))
      reflection_list = reflection_list.select(perm)

    for refl in reflection_list:
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

  try:
    import matplotlib

    if not params.output.show_plot:
      # http://matplotlib.org/faq/howto_faq.html#generate-images-without-having-a-window-appear
      matplotlib.use('Agg') # use a non-interactive backend
    from matplotlib import pyplot
  except ImportError:
    raise Sorry("matplotlib must be installed to generate a plot.")

  fig = pyplot.figure()
  fig.set_size_inches(params.output.size_inches)
  fig.set_dpi(params.output.dpi)
  pyplot.axes().set_aspect('equal')
  pyplot.scatter(obs_x, obs_y, marker='x', c='red', s=10, alpha=1)
  pyplot.scatter(pred_x, pred_y, marker='+', c='blue')
  assert len(detector) == 1
  panel = detector[0]
  pyplot.xlim(0, panel.get_image_size_mm()[0])
  pyplot.ylim(0, panel.get_image_size_mm()[1])
  pyplot.gca().invert_yaxis()
  pyplot.title('Centroid x,y-coordinates')
  pyplot.xlabel('x-coordinate (mm)')
  pyplot.ylabel('y-coordinate (mm)')
  if params.output.file_name is not None:
    pyplot.savefig(params.output.file_name,
                   size_inches=params.output.size_inches,
                   dpi=params.output.dpi,
                   bbox_inches='tight')
  if params.output.show_plot:
    pyplot.show()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
