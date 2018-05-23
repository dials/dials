# LIBTBX_SET_DISPATCHER_NAME dials.show

from __future__ import absolute_import, division, print_function

import iotbx.phil
import numpy
from libtbx import table_utils

help_message = '''

Examples::

  dials.show datablock.json

  dials.show experiments.json

  dials.show image_*.cbf

  dials.show reflections.pickle

'''

phil_scope = iotbx.phil.parse("""\
show_scan_varying = False
  .type = bool
  .help = "Whether or not to show the crystal at each scan point."
show_all_reflection_data = False
  .type = bool
  .help = "Whether or not to print individual reflections"
show_intensities = False
  .type = bool
show_centroids = False
  .type = bool
show_profile_fit = False
  .type = bool
show_flags = False
  .type = bool
  .help = "Show a summary table of reflection flags"
max_reflections = None
  .type = int
  .help = "Limit the number of reflections in the output."
show_panel_distance = False
  .type = bool
  .help = "Show distance to individual panels along normal."
""", process_includes=True)


def beam_centre_mm(detector, s0):
  x, y = (None, None)
  for panel_id, panel in enumerate(detector):
    try:
      x, y = panel.get_ray_intersection(s0)
    except RuntimeError:
      continue
    else:
      if panel.is_coord_valid_mm((x, y)):
        break
      else:
        x, y = (None, None)
  return panel_id, (x, y)


def show_beam(detector, beam):

  # standard static beam model string
  s = str(beam)

  # report whether the beam is scan-varying
  if beam.num_scan_points > 0:
    s += "    s0 sampled at " + str(beam.num_scan_points) + " scan points\n"

  # add static model beam centres
  panel_id, (x, y) = beam_centre_mm(detector, beam.get_s0())
  if panel_id >= 0 and x is not None and y is not None:
    x_px, y_px = detector[panel_id].millimeter_to_pixel((x, y))
    if len(detector) > 1:
      beam_centre_mm_str = "Beam centre (mm): panel %i, (%.2f,%.2f)" %(
        panel_id, x, y)
      beam_centre_px_str = "Beam centre (px): panel %i, (%.2f,%.2f)" %(
        panel_id, x_px, y_px)
    else:
      beam_centre_mm_str = "Beam centre (mm): (%.2f,%.2f)" %(x, y)
      beam_centre_px_str = "Beam centre (px): (%.2f,%.2f)" %(x_px, y_px)

    s += beam_centre_mm_str + '\n' + beam_centre_px_str + '\n'

  # report range of scan-varying model beam centres
  if beam.num_scan_points > 0:
    # get scan-varying beam centres, ensuring all on same panel
    sv_s0 = beam.get_s0_at_scan_points()
    impacts = [beam_centre_mm(detector, s0) for s0 in sv_s0]
    pnl, xy = zip(*impacts)
    uniq_pnls = set(pnl)
    if len(uniq_pnls) > 1 or min(uniq_pnls) < 0: return static_str + sv_str
    pnl = list(uniq_pnls)[0]
    x_mm, y_mm = zip(*xy)

    # convert to pixels
    xy = [detector[pnl].millimeter_to_pixel(e) for e in xy]
    x_px, y_px = zip(*xy)

    s += "Beam centre range (mm): ([%.2f,%.2f],[%.2f,%.2f])\n" % (min(x_mm),
        max(x_mm), min(y_mm), max(y_mm))
    s += "Beam centre range (px): ([%.2f,%.2f],[%.2f,%.2f])\n" % (min(x_px),
        max(x_px), min(y_px), max(y_px))

  return s

def show_goniometer(goniometer):

  # standard static goniometer model string
  s = str(goniometer)

  # report whether the goniometer is scan-varying
  if goniometer.num_scan_points > 0:
     s += ("    Setting rotation sampled at " +
     str(goniometer.num_scan_points) + " scan points\n")

  return s

def run(args):
  import dials.util.banner
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] datablock.json | experiments.json | image_*.cbf" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_datablocks=True,
    read_datablocks_from_images=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  experiments = flatten_experiments(params.input.experiments)
  datablocks = flatten_datablocks(params.input.datablock)
  reflections = flatten_reflections(params.input.reflections)

  if len(datablocks) == 0 and len(experiments) == 0 and len(reflections) == 0:
    parser.print_help()
    exit()

  if len(experiments):
    print(show_experiments(
      experiments, show_panel_distance=params.show_panel_distance,
      show_scan_varying=params.show_scan_varying))

  if len(datablocks):
    print(show_datablocks(
      datablocks, show_panel_distance=params.show_panel_distance))

  if len(reflections):
    print(show_reflections(
      reflections, show_intensities=params.show_intensities,
      show_profile_fit=params.show_profile_fit,
      show_centroids=params.show_centroids,
      show_all_reflection_data=params.show_all_reflection_data,
      show_flags=params.show_flags,
      max_reflections=params.max_reflections))


def show_experiments(experiments, show_panel_distance=False,
                     show_scan_varying=False):

  text = []

  for i_expt, expt in enumerate(experiments):
    text.append("Experiment %i:" %i_expt)
    text.append(str(expt.detector))
    text.append('Max resolution (at corners): %f' % (
      expt.detector.get_max_resolution(expt.beam.get_s0())))
    text.append('Max resolution (inscribed):  %f' % (
      expt.detector.get_max_inscribed_resolution(expt.beam.get_s0())))
    if show_panel_distance:
      for ipanel, panel in enumerate(expt.detector):
        from scitbx import matrix
        fast = matrix.col(panel.get_fast_axis())
        slow = matrix.col(panel.get_slow_axis())
        normal = fast.cross(slow).normalize()
        origin = matrix.col(panel.get_origin())
        distance = origin.dot(normal)
        fast_origin = - (origin - distance * normal).dot(fast)
        slow_origin = - (origin - distance * normal).dot(slow)
        text.append('Panel %d: distance %.2f origin %.2f %.2f' % \
          (ipanel, distance, fast_origin, slow_origin))
      text.append('')
    text.append('')
    text.append(show_beam(expt.detector, expt.beam))
    if expt.scan is not None:
      text.append(str(expt.scan))
    if expt.goniometer is not None:
      text.append(show_goniometer(expt.goniometer))
    from cStringIO import StringIO
    s = StringIO()
    expt.crystal.show(show_scan_varying=show_scan_varying, out=s)
    text.append(s.getvalue())
    if expt.crystal.num_scan_points:
      from scitbx.array_family import flex
      from cctbx import uctbx
      abc = flex.vec3_double()
      angles = flex.vec3_double()
      for n in range(expt.crystal.num_scan_points):
        a, b, c, alpha, beta, gamma = expt.crystal.get_unit_cell_at_scan_point(n).parameters()
        abc.append((a, b, c))
        angles.append((alpha, beta, gamma))
      a, b, c = abc.mean()
      alpha, beta, gamma = angles.mean()
      mean_unit_cell = uctbx.unit_cell((a, b, c, alpha, beta, gamma))
      text.append('  Average unit cell: %s' %mean_unit_cell)
    #text.append('')
    if expt.profile is not None:
      text.append(str(expt.profile))
    if expt.scaling_model is not None:
      text.append(str(expt.scaling_model.show()))
  return '\n'.join(text)


def show_datablocks(datablocks, show_panel_distance=False):
  text = []
  for datablock in datablocks:
    if datablock.format_class() is not None:
      text.append('Format: %s' %datablock.format_class())
    imagesets = datablock.extract_imagesets()
    for imageset in imagesets:
      try: text.append(imageset.get_template())
      except Exception: pass
      detector = imageset.get_detector()
      if detector is not None:
        text.append(str(detector))
        text.append('Max resolution (at corners): %f' % (
          detector.get_max_resolution(imageset.get_beam().get_s0())))
        text.append('Max resolution (inscribed):  %f' % (
          detector.get_max_inscribed_resolution(imageset.get_beam().get_s0())))
        if show_panel_distance:
          for ipanel, panel in enumerate(detector):
            from scitbx import matrix
            fast = matrix.col(panel.get_fast_axis())
            slow = matrix.col(panel.get_slow_axis())
            normal = fast.cross(slow)
            origin = matrix.col(panel.get_origin())
            distance = origin.dot(normal)
            fast_origin = - (origin - distance * normal).dot(fast)
            slow_origin = - (origin - distance * normal).dot(slow)
            text.append('Panel %d: distance %.2f origin %.2f %.2f' % \
              (ipanel, distance, fast_origin, slow_origin))
          text.append('')
        text.append('')
        text.append(show_beam(detector, imageset.get_beam()))
      if imageset.get_scan() is not None:
        text.append(str(imageset.get_scan()))
      if imageset.get_goniometer() is not None:
        text.append(str(imageset.get_goniometer()))
  return '\n'.join(text)


def _create_flag_count_table(table):
  """Generate a summary table of flag values in a reflection table.

  :param table: A reflection table
  :returns:     A string of the formatted flags table
  """

  # Calculate the counts of entries that match each flag
  numpy_flags = table["flags"].as_numpy_array()
  flag_count = {flag: numpy.sum(numpy_flags & value != 0) for value, flag in table.flags.values.items()}

  # Work out the numeric-value order of the flags
  flag_order = sorted(table.flags.values.values(), key=lambda x: x.real)

  # Build the actual table
  flag_rows = [["Flag", "Count", "%"]]
  max_count_len = max(5, len(str(max(flag_count.values()))))
  last_flag = None
  for flag in flag_order:
    indent = ""
    # As a hint for reading, indent any 'summary' flags.
    # A summary flag is any flag which overlaps with the previous one.
    if last_flag and (last_flag.real & flag.real):
      indent = "  "
    last_flag = flag
    # Add the row to the table we're building
    flag_rows.append([indent + flag.name,
                      "{:{:d}d}".format(flag_count[flag], max_count_len),
                      "{:5.01f}".format(100*flag_count[flag] / len(table))])

  # Build the array of output strings
  text = []
  text.append("Reflection flags:")
  text.append(
    table_utils.format(flag_rows, has_header=True, prefix="| ", postfix=" |"))
  return "\n".join(text)


def show_reflections(reflections, show_intensities=False, show_profile_fit=False,
                     show_centroids=False, show_all_reflection_data=False,
                     show_flags=False, max_reflections=None):

  text = []

  import collections
  from libtbx.containers import OrderedSet
  formats = collections.OrderedDict((
    ('miller_index', '%i, %i, %i'),
    ('d','%.2f'),
    ('qe','%.3f'),
    ('dqe','%.3f'),
    ('id','%i'),
    ('imageset_id','%i'),
    ('panel','%i'),
    ('flags', '%i'),
    ('background.mean', '%.1f'),
    ('background.dispersion','%.1f'),
    ('background.mse', '%.1f'),
    ('background.sum.value', '%.1f'),
    ('background.sum.variance', '%.1f'),
    ('intensity.prf.value','%.1f'),
    ('intensity.prf.variance','%.1f'),
    ('intensity.sum.value','%.1f'),
    ('intensity.sum.variance','%.1f'),
    ('intensity.cor.value','%.1f'),
    ('intensity.cor.variance','%.1f'),
    ('lp','%.3f'),
    ('num_pixels.background','%i'),
    ('num_pixels.background_used','%i'),
    ('num_pixels.foreground','%i'),
    ('num_pixels.valid','%i'),
    ('partial_id','%i'),
    ('partiality','%.4f'),
    ('profile.correlation','%.3f'),
    ('profile.rmsd','%.3f'),
    ('xyzcal.mm','%.2f, %.2f, %.2f'),
    ('xyzcal.px','%.2f, %.2f, %.2f'),
    ('delpsical.rad','%.3f'),
    ('delpsical2','%.3f'),
    ('delpsical.weights','%.3f'),
    ('xyzobs.mm.value','%.2f, %.2f, %.2f'),
    ('xyzobs.mm.variance','%.4e, %.4e, %.4e'),
    ('xyzobs.px.value','%.2f, %.2f, %.2f'),
    ('xyzobs.px.variance','%.4f, %.4f, %.4f'),
    ('s1','%.4f, %.4f, %.4f'),
    ('shoebox','%.1f'),
    ('rlp','%.4f, %.4f, %.4f'),
    ('zeta','%.3f'),
    ('x_resid','%.3f'),
    ('x_resid2','%.3f'),
    ('y_resid','%.3f'),
    ('y_resid2','%.3f'),
    ('kapton_absorption_correction','%.3f'),
    ('kapton_absorption_correction_sigmas','%.3f'),
    ))

  for rlist in reflections:
    from dials.array_family import flex
    from dials.algorithms.shoebox import MaskCode
    foreground_valid = MaskCode.Valid | MaskCode.Foreground
    text.append('')
    text.append('Reflection list contains %i reflections' %(len(rlist)))

    if len(rlist) == 0:
      continue

    rows = [["Column", "min", "max", "mean"]]
    for k, col in rlist.cols():
      if k in formats and not "%" in formats[k]:
        # Allow blanking out of entries that wouldn't make sense
        rows.append([k, formats[k], formats[k], formats[k]])
      elif type(col) in (flex.double, flex.int, flex.size_t):
        if type(col) in (flex.int, flex.size_t):
          col = col.as_double()
        rows.append([k, formats[k] %flex.min(col), formats[k] %flex.max(col),
                     formats[k]%flex.mean(col)])
      elif type(col) in (flex.vec3_double, flex.miller_index):
        if isinstance(col, flex.miller_index):
          col = col.as_vec3_double()
        rows.append([k, formats[k] %col.min(), formats[k] %col.max(),
                     formats[k]%col.mean()])
      elif isinstance(col, flex.shoebox):
        rows.append([k, "", "", ""])
        si = col.summed_intensity().observed_value()
        rows.append(["  summed I", formats[k] %flex.min(si), formats[k] %flex.max(si),
                     formats[k]%flex.mean(si)])
        x1, x2, y1, y2, z1, z2 = col.bounding_boxes().parts()
        bbox_sizes = ((z2-z1)*(y2-y1)*(x2-x1)).as_double()
        rows.append(["  N pix", formats[k] %flex.min(bbox_sizes), formats[k] %flex.max(bbox_sizes),
                     formats[k]%flex.mean(bbox_sizes)])
        fore_valid = col.count_mask_values(foreground_valid).as_double()
        rows.append(["  N valid foreground pix", formats[k] %flex.min(fore_valid), formats[k] %flex.max(fore_valid),
                     formats[k]%flex.mean(fore_valid)])

    text.append(
      table_utils.format(rows, has_header=True, prefix="| ", postfix=" |"))

    if show_flags:
      text.append(_create_flag_count_table(rlist))

  intensity_keys = (
    'miller_index', 'd', 'intensity.prf.value', 'intensity.prf.variance',
    'intensity.sum.value', 'intensity.sum.variance', 'background.mean',
    'profile.correlation', 'profile.rmsd'
  )

  profile_fit_keys = ('miller_index', 'd',)

  centroid_keys = (
    'miller_index', 'd', 'xyzcal.mm', 'xyzcal.px', 'xyzobs.mm.value',
    'xyzobs.mm.variance', 'xyzobs.px.value', 'xyzobs.px.variance'
  )

  keys_to_print = OrderedSet()

  if show_intensities:
    for k in intensity_keys: keys_to_print.add(k)
  if show_profile_fit:
    for k in profile_fit_keys: keys_to_print.add(k)
  if show_centroids:
    for k in centroid_keys: keys_to_print.add(k)
  if show_all_reflection_data:
    for k in formats: keys_to_print.add(k)

  def format_column(key, data, format_strings=None):
    if isinstance(data, flex.vec3_double):
      c_strings = [c.as_string(format_strings[i].strip()) for i, c in enumerate(data.parts())]
    elif isinstance(data, flex.miller_index):
      c_strings = [c.as_string(format_strings[i].strip()) for i, c in enumerate(data.as_vec3_double().parts())]
    elif isinstance(data, flex.size_t):
      c_strings = [data.as_int().as_string(format_strings[0].strip())]
    elif isinstance(data, flex.shoebox):
      x1, x2, y1, y2, z1, z2 = data.bounding_boxes().parts()
      bbox_sizes = ((z2-z1)*(y2-y1)*(x2-x1)).as_double()
      c_strings = [bbox_sizes.as_string(format_strings[0].strip())]
      key += " (N pix)"
    else:
      c_strings = [data.as_string(format_strings[0].strip())]

    column = flex.std_string()
    max_element_lengths = [c.max_element_length() for c in c_strings]
    for i in range(len(c_strings[0])):

      column.append(('%%%is' %len(key)) %', '.join(
        ('%%%is' %max_element_lengths[j]) %c_strings[j][i]
        for j in range(len(c_strings))))
    return column


  if keys_to_print:
    keys = [k for k in keys_to_print if k in rlist]
    rows = [keys]
    if max_reflections is not None:
      max_reflections = min(len(rlist), max_reflections)

    columns = []

    for k in keys:
      columns.append(format_column(k, rlist[k], format_strings=formats[k].split(',')))

    text.append('')
    text.append('Printing %i of %i reflections:' %(max_reflections, len(rlist)))
    line = []
    for j in range(len(columns)):
      key = keys[j]
      if key == 'shoebox': key += " (N pix)"
      width = max(len(key), columns[j].max_element_length())
      line.append('%%%is' %width %key)
    text.append(' '.join(line))
    for i in range(max_reflections):
      line = []
      for j in range(len(columns)):
        line.append(columns[j][i])
      text.append(' '.join(line))

  return '\n'.join(text)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
