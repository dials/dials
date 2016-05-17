# LIBTBX_SET_DISPATCHER_NAME dials.show

from __future__ import division

import iotbx.phil

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
max_reflections = None
  .type = int
  .help = "Limit the number of reflections in the output."
""", process_includes=True)


def beam_centre(detector, beam):
  s0 = beam.get_s0()
  x, y = (None, None)
  for panel_id, panel in enumerate(detector):
    try:
      x, y = panel.get_ray_intersection(s0)
    except RuntimeError:
      continue
    else:
      if panel.is_coord_valid_mm((x, y)):
        break
  return panel_id, (x, y)


def run(args):

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

  for i_expt, expt in enumerate(experiments):
    print "Experiment %i:" %i_expt
    print str(expt.detector)
    print 'Max resolution (at corners): %f' % (
      expt.detector.get_max_resolution(expt.beam.get_s0()))
    print 'Max resolution (inscribed):  %f' % (
      expt.detector.get_max_inscribed_resolution(expt.beam.get_s0()))
    print ''
    panel_id, (x, y) = beam_centre(expt.detector, expt.beam)
    if panel_id >= 0 and x is not None and y is not None:
      if len(expt.detector) > 1:
        beam_centre_str = "Beam centre: panel %i, (%.2f,%.2f)" %(panel_id, x, y)
      else:
        beam_centre_str = "Beam centre: (%.2f,%.2f)" %(x, y)
    else:
      beam_centre_str = ""
    print str(expt.beam) + beam_centre_str + '\n'
    if expt.scan is not None:
      print expt.scan
    if expt.goniometer is not None:
      print expt.goniometer
    expt.crystal.show(show_scan_varying=params.show_scan_varying)
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
      print "  Average unit cell: %s" %mean_unit_cell
    print

  for datablock in datablocks:
    if datablock.format_class() is not None:
      print 'Format: %s' %datablock.format_class()
    imagesets = datablock.extract_imagesets()
    for imageset in imagesets:
      try: print imageset.get_template()
      except Exception: pass
      detector = imageset.get_detector()
      print str(detector) + 'Max resolution: %f\n' %(
        detector.get_max_resolution(imageset.get_beam().get_s0()))
      panel_id, (x, y) = beam_centre(detector, imageset.get_beam())
      if panel_id >= 0 and x is not None and y is not None:
        if len(detector) > 1:
          beam_centre_str = "Beam centre: panel %i, (%.2f,%.2f)" %(panel_id, x, y)
        else:
          beam_centre_str = "Beam centre: (%.2f,%.2f)" %(x, y)
      else:
        beam_centre_str = ""
      print str(imageset.get_beam()) + beam_centre_str + '\n'
      if imageset.get_scan() is not None:
        print imageset.get_scan()
      if imageset.get_goniometer() is not None:
        print imageset.get_goniometer()

  from libtbx.containers import OrderedDict, OrderedSet
  formats = OrderedDict([
    ('miller_index', '%i, %i, %i'),
    ('d','%.2f'),
    ('dqe','%.3f'),
    ('id','%i'),
    ('imageset_id','%i'),
    ('panel','%i'),
    ('flags', '%i'),
    ('background.mean', '%.1f'),
    ('intensity.prf.value','%.1f'),
    ('intensity.prf.variance','%.1f'),
    ('intensity.sum.value','%.1f'),
    ('intensity.sum.variance','%.1f'),
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
    ('xyzobs.mm.value','%.2f, %.2f, %.2f'),
    ('xyzobs.mm.variance','%.4e, %.4e, %.4e'),
    ('xyzobs.px.value','%.2f, %.2f, %.2f'),
    ('xyzobs.px.variance','%.4f, %.4f, %.4f'),
    ('s1','%.4f, %.4f, %.4f'),
    ('rlp','%.4f, %.4f, %.4f'),
    ('zeta','%.3f'),
    ('x_resid','%.3f'),
    ('x_resid2','%.3f'),
    ('y_resid','%.3f'),
    ('y_resid2','%.3f'),
    ])

  for rlist in reflections:
    from cctbx.array_family import flex
    print
    print "Reflection list contains %i reflections" %(len(rlist))
    rows = [["Column", "min", "max", "mean"]]
    for k, col in rlist.cols():
      if type(col) in (flex.double, flex.int, flex.size_t):
        if type(col) in (flex.int, flex.size_t):
          col = col.as_double()
        rows.append([k, formats[k] %flex.min(col), formats[k] %flex.max(col),
                     formats[k]%flex.mean(col)])
      elif type(col) in (flex.vec3_double, flex.miller_index):
        if type(col) == flex.miller_index:
          col = col.as_vec3_double()
        rows.append([k, formats[k] %col.min(), formats[k] %col.max(),
                     formats[k]%col.mean()])

    from libtbx import table_utils
    print table_utils.format(rows, has_header=True, prefix="| ", postfix=" |")

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

  if params.show_intensities:
    for k in intensity_keys: keys_to_print.add(k)
  if params.show_profile_fit:
    for k in profile_fit_keys: keys_to_print.add(k)
  if params.show_centroids:
    for k in centroid_keys: keys_to_print.add(k)
  if params.show_all_reflection_data:
    for k in formats: keys_to_print.add(k)

  def format_column(key, data, format_strings=None):
    if isinstance(data, flex.vec3_double):
      c_strings = [c.as_string(format_strings[i].strip()) for i, c in enumerate(data.parts())]
    elif isinstance(data, flex.miller_index):
      c_strings = [c.as_string(format_strings[i].strip()) for i, c in enumerate(data.as_vec3_double().parts())]
    elif isinstance(data, flex.size_t):
      c_strings = [data.as_int().as_string(format_strings[0].strip())]
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
    max_reflections = len(rlist)
    if params.max_reflections is not None:
      max_reflections = min(len(rlist), params.max_reflections)

    columns = []

    for k in keys:
      columns.append(format_column(k, rlist[k], format_strings=formats[k].split(',')))

    print
    print "Printing %i of %i reflections:" %(max_reflections, len(rlist))
    for j in range(len(columns)):
      key = keys[j]
      width = max(len(key), columns[j].max_element_length())
      print ("%%%is" %width) %key,
    print
    for i in range(max_reflections):
      for j in range(len(columns)):
        print columns[j][i],
      print


    #for i in range(max_reflections):
      #r = rlist[i]
      #rows.append([formats[k] %r[k] for k in keys])

    #print "Printing %i of %i reflections:" %(max_reflections, len(rlist))
    ##print table_utils.format(rows, has_header=True, prefix="| ", postfix=" |")
    ##print table_utils.format(rows, has_header=False, prefix="", postfix="",
                             ##delim=" ")
    ##from dials.util.tabulate import tabulate
    ##print tabulate(rows)
    #import asciitable
    #asciitable.write(rows)

  return

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
