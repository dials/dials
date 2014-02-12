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
  from scitbx import matrix
  from libtbx.phil import command_line
  importer = Importer(args, check_format=False)
  assert len(importer.reflections) == 1
  assert len(importer.datablocks) == 1
  reflections = importer.reflections[0]
  imageset = importer.datablocks[0].extract_imagesets()[0]
  detector = imageset.get_detector()
  scan = imageset.get_scan()
  args = importer.unhandled_arguments
  cmd_line = command_line.argument_interpreter(master_params=master_phil_scope)
  working_phil = cmd_line.process_and_fetch(args=args)
  params = working_phil.extract()
  observed_xy = flex.vec2_double()
  predicted_xy = flex.vec2_double()
  panel_origin_shifts = {0: (0,0,0)}
  try:
    hierarchy = detector.hierarchy()
  except AttributeError, e:
    hierarchy = None
  for i_panel in range(1, len(detector)):
    origin_shift = matrix.col(detector[0].get_origin()) \
      - matrix.col(detector[i_panel].get_origin())
    panel_origin_shifts[i_panel] = origin_shift
  from dials.model.data import ReflectionList
  reflection_list = ReflectionList.from_table(reflections)
  if len(params.scan_range):
    sel = flex.bool(reflection_list.size(), False)
    centroid_positions = reflection_list.centroid_position()
    if (centroid_positions.norms().all_eq(0) or
        not reflection_list.miller_index().all_eq((0,0,0))):
      centroids_frame = reflection_list.frame_number()
    else:
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
      if refl.miller_index is not None and refl.image_coord_mm != (0,0):
        x, y = refl.image_coord_mm
      else:
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
        if refl.panel_number > 0:
          s1 = detector[refl.panel_number].get_lab_coord(centroid_position[:2])
          if hierarchy is not None:
            centroid_position = hierarchy.get_ray_intersection(s1)
          else:
            centroid_position = detector[0].get_ray_intersection(s1)
        x, y = centroid_position[:2]
      observed_xy.append((x,y))
    elif refl.image_coord_px != (0, 0):
      x, y = refl.image_coord_mm
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
  pyplot.scatter(obs_x, obs_y, marker='o', c='white', s=10, alpha=1)
  pyplot.scatter(pred_x, pred_y, marker='+', c='blue')
  #assert len(detector) == 1
  panel = detector[0]
  #if len(detector) > 1:
  xmin = max([detector[i_panel].get_image_size_mm()[0] + panel_origin_shifts[i_panel][0]
              for i_panel in range(len(detector))])
  xmax = max([detector[i_panel].get_image_size_mm()[0] + panel_origin_shifts[i_panel][0]
              for i_panel in range(len(detector))])
  ymax = max([detector[i_panel].get_image_size_mm()[1] + panel_origin_shifts[i_panel][1]
              for i_panel in range(len(detector))])
  ymax = max([detector[i_panel].get_image_size_mm()[1] + panel_origin_shifts[i_panel][1]
              for i_panel in range(len(detector))])
  if hierarchy is not None:
    beam_centre = hierarchy.get_beam_centre(imageset.get_beam().get_s0())
  else:
    beam_centre = detector[0].get_beam_centre(imageset.get_beam().get_s0())
  pyplot.scatter([beam_centre[0]], [beam_centre[1]], marker='+', c='blue', s=100)
  pyplot.xlim(0, xmax)
  pyplot.ylim(0, ymax)
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
