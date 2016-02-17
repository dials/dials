# LIBTBX_SET_DISPATCHER_NAME dev.dials.plot_find_spots_client
from __future__ import division

from cctbx.array_family import flex
import iotbx.phil

help_message = '''
'''

phil_scope = iotbx.phil.parse("""\
grid = None
  .type = ints(size=2, value_min=1)
stereographic_projections = False
  .type = bool
positions = None
  .type = path
cmap = YlOrRd
  .type = str
invalid = white
  .type = str
""", process_includes=True)


def run(args):
  from dials.util.options import OptionParser
  import libtbx.load_env

  usage = "%s [options] find_spots.json" %(
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    epilog=help_message)

  params, options, args = parser.parse_args(
    show_diff_phil=True, return_unhandled=True)

  positions = None
  if params.positions is not None:
    with open(params.positions, 'rb') as f:
      positions = flex.vec2_double()
      for line in f.readlines():
        line = line.replace('(', ' ').replace(')', '').replace(',', ' ').strip().split()
        assert len(line) == 3
        i, x, y = [float(l) for l in line]
        positions.append((x, y))

  assert len(args) == 1
  json_file = args[0]
  import json

  with open(json_file, 'rb') as f:
    results = json.load(f)

  n_indexed = flex.double()
  fraction_indexed = flex.double()
  n_spots = flex.double()
  n_lattices = flex.double()
  crystals = []
  image_names = flex.std_string()

  for r in results:
    n_spots.append(r['n_spots_total'])
    image_names.append(str(r['image']))
    if 'n_indexed' in r:
      n_indexed.append(r['n_indexed'])
      n_lattices.append(len(r['lattices']))
      for d in r['lattices']:
        from dxtbx.serialize.crystal import from_dict
        crystals.append(from_dict(d['crystal']))
    else:
      n_indexed.append(0)
      n_lattices.append(0)

  if n_indexed.size():
    sel = n_spots > 0
    fraction_indexed = flex.double(n_indexed.size(), 0)
    fraction_indexed.set_selected(
      sel, n_indexed.select(sel)/n_spots.select(sel))

  import matplotlib
  matplotlib.use('Agg')
  from matplotlib import pyplot

  blue = '#3498db'
  red = '#e74c3c'

  marker = 'o'
  alpha = 0.5
  lw = 0

  plot = True
  table = True
  grid = params.grid

  from libtbx import group_args
  from dials.algorithms.peak_finding.per_image_analysis \
       import plot_stats, print_table

  estimated_d_min = flex.double()
  d_min_distl_method_1 = flex.double()
  d_min_distl_method_2 = flex.double()
  n_spots_total = flex.int()
  n_spots_no_ice = flex.int()
  total_intensity = flex.double()

  for d in results:
    estimated_d_min.append(d['estimated_d_min'])
    d_min_distl_method_1.append(d['d_min_distl_method_1'])
    d_min_distl_method_2.append(d['d_min_distl_method_2'])
    n_spots_total.append(d['n_spots_total'])
    n_spots_no_ice.append(d['n_spots_no_ice'])
    total_intensity.append(d['total_intensity'])

  stats = group_args(image=image_names,
                     n_spots_total=n_spots_total,
                     n_spots_no_ice=n_spots_no_ice,
                     n_spots_4A=None,
                     n_indexed=n_indexed,
                     fraction_indexed=fraction_indexed,
                     total_intensity=total_intensity,
                     estimated_d_min=estimated_d_min,
                     d_min_distl_method_1=d_min_distl_method_1,
                     d_min_distl_method_2=d_min_distl_method_2,
                     noisiness_method_1=None,
                     noisiness_method_2=None)

  if plot:
    plot_stats(stats)
    pyplot.clf()
  if table:
    print_table(stats)

  print "Number of indexed lattices: ", (n_indexed > 0).count(True)

  print "Number with valid d_min but failed indexing: ", (
    (d_min_distl_method_1 > 0) &
    (d_min_distl_method_2 > 0) &
    (estimated_d_min > 0) &
    (n_indexed == 0)).count(True)

  n_rows = 10
  n_rows = min(n_rows, len(n_spots_total))
  perm_n_spots_total = flex.sort_permutation(n_spots_total, reverse=True)
  print 'Top %i images sorted by number of spots:' %n_rows
  print_table(stats, perm=perm_n_spots_total, n_rows=n_rows)

  n_bins = 20
  spot_count_histogram(
    n_spots_total, n_bins=n_bins, filename='hist_n_spots_total.png', log=True)
  spot_count_histogram(
    n_spots_no_ice, n_bins=n_bins, filename='hist_n_spots_no_ice.png', log=True)
  spot_count_histogram(
    n_indexed.select(n_indexed > 0), n_bins=n_bins, filename='hist_n_indexed.png', log=False)

  if len(crystals):
    plot_unit_cell_histograms(crystals)

  if params.stereographic_projections and len(crystals):
    from dxtbx.datablock import DataBlockFactory
    datablocks = DataBlockFactory.from_filenames(
      [image_names[0]], verbose=False)
    assert len(datablocks) == 1
    imageset = datablocks[0].extract_imagesets()[0]
    s0 = imageset.get_beam().get_s0()
    # XXX what if no goniometer?
    rotation_axis = imageset.get_goniometer().get_rotation_axis()

    indices = ((1,0,0), (0,1,0), (0,0,1))
    for i, index in enumerate(indices):

      from cctbx import crystal, miller
      from scitbx import matrix
      miller_indices = flex.miller_index([index])
      symmetry = crystal.symmetry(
        unit_cell=crystals[0].get_unit_cell(),
        space_group=crystals[0].get_space_group())
      miller_set = miller.set(symmetry, miller_indices)
      d_spacings = miller_set.d_spacings()
      d_spacings = d_spacings.as_non_anomalous_array().expand_to_p1()
      d_spacings = d_spacings.generate_bijvoet_mates()
      miller_indices = d_spacings.indices()

      # plane normal
      d0 = matrix.col(s0).normalize()
      d1 = d0.cross(matrix.col(rotation_axis)).normalize()
      d2 = d1.cross(d0).normalize()
      reference_poles = (d0, d1, d2)

      from dials.command_line.stereographic_projection import stereographic_projection
      projections = []

      for cryst in crystals:
        reciprocal_space_points = list(cryst.get_U() * cryst.get_B()) * miller_indices.as_vec3_double()
        projections.append(stereographic_projection(
          reciprocal_space_points, reference_poles))

        #from dials.algorithms.indexing.compare_orientation_matrices import \
        #  difference_rotation_matrix_and_euler_angles
        #R_ij, euler_angles, cb_op = difference_rotation_matrix_and_euler_angles(
        #  crystals[0], cryst)
        #print max(euler_angles)

      from dials.command_line.stereographic_projection import plot_projections
      plot_projections(projections, filename='projections_%s.png' %('hkl'[i]))
      pyplot.clf()

  def plot_grid(values, grid, file_name, cmap=pyplot.cm.Reds,
                vmin=None, vmax=None, invalid='white'):
    values = values.as_double()
    # At DLS, fast direction appears to be largest direction
    if grid[0] > grid[1]:
      values.reshape(flex.grid(reversed(grid)))
      values = values.matrix_transpose()
    else:
      values.reshape(flex.grid(grid))

    Z = values.as_numpy_array()

    #f, (ax1, ax2) = pyplot.subplots(2)
    f, ax1 = pyplot.subplots(1)

    mesh1 = ax1.pcolormesh(
      values.as_numpy_array(), cmap=cmap, vmin=vmin, vmax=vmax)
    mesh1.cmap.set_under(color=invalid, alpha=None)
    mesh1.cmap.set_over(color=invalid, alpha=None)
    #mesh2 = ax2.contour(Z, cmap=cmap, vmin=vmin, vmax=vmax)
    #mesh2 = ax2.contourf(Z, cmap=cmap, vmin=vmin, vmax=vmax)
    ax1.set_aspect('equal')
    ax1.invert_yaxis()
    #ax2.set_aspect('equal')
    #ax2.invert_yaxis()
    pyplot.colorbar(mesh1, ax=ax1)
    #pyplot.colorbar(mesh2, ax=ax2)
    pyplot.savefig(file_name, dpi=600)
    pyplot.clf()

  def plot_positions(values, positions, file_name, cmap=pyplot.cm.Reds,
                     vmin=None, vmax=None, invalid='white'):
    values = values.as_double()
    assert positions.size() >= values.size()
    positions = positions[:values.size()]

    if vmin is None:
      vmin = flex.min(values)
    if vmax is None:
      vmax = flex.max(values)

    x, y = positions.parts()
    dx = flex.abs(x[1:] - x[:-1])
    dy = flex.abs(y[1:] - y[:-1])
    dx = dx.select(dx > 0)
    dy = dy.select(dy > 0)

    scale = 1/flex.min(dx)
    #print scale
    x = (x * scale).iround()
    y = (y * scale).iround()

    from libtbx.math_utils import iceil
    z = flex.double(flex.grid(iceil(flex.max(y))+1, iceil(flex.max(x))+1), -2)
    #print z.all()
    for x_, y_, z_ in zip(x, y, values):
      z[y_, x_] = z_

    plot_grid(z.as_1d(), z.all(), file_name, cmap=cmap, vmin=vmin, vmax=vmax,
              invalid=invalid)
    return

  if grid is not None or positions is not None:
    if grid is not None:
      positions = tuple(reversed(grid))
      plotter = plot_grid
    else:
      plotter = plot_positions

    cmap = pyplot.get_cmap(params.cmap)
    plotter(n_spots_total, positions, 'grid_spot_count_total.png', cmap=cmap,
            invalid=params.invalid)
    plotter(n_spots_no_ice, positions, 'grid_spot_count_no_ice.png', cmap=cmap,
            invalid=params.invalid)
    plotter(total_intensity, positions, 'grid_total_intensity.png', cmap=cmap,
            invalid=params.invalid)
    if flex.max(n_indexed) > 0:
      plotter(n_indexed, positions, 'grid_n_indexed.png', cmap=cmap,
              invalid=params.invalid)
      plotter(fraction_indexed, positions, 'grid_fraction_indexed.png',
              cmap=cmap, vmin=0, vmax=1, invalid=params.invalid)

    for i, d_min in enumerate((estimated_d_min, d_min_distl_method_1, d_min_distl_method_2)):
      from cctbx import uctbx
      d_star_sq = uctbx.d_as_d_star_sq(d_min)
      d_star_sq.set_selected(d_star_sq == 1, 0)
      vmin = flex.min(d_star_sq.select(d_star_sq > 0))
      vmax = flex.max(d_star_sq)

      vmin = flex.min(d_min.select(d_min > 0))
      vmax = flex.max(d_min)
      cmap = pyplot.get_cmap('%s_r' %params.cmap)
      d_min.set_selected(d_min <= 0, vmax)

      if i == 0:
        plotter(d_min, positions, 'grid_d_min.png', cmap=cmap, vmin=vmin,
                vmax=vmax, invalid=params.invalid)
      else:
        plotter(
          d_min, positions, 'grid_d_min_method_%i.png' %i, cmap=cmap,
          vmin=vmin, vmax=vmax, invalid=params.invalid)

  if flex.max(n_indexed) > 0:
    pyplot.hexbin(
      n_spots, n_indexed, bins='log', cmap=pyplot.cm.jet, gridsize=50)
    pyplot.colorbar()
    #pyplot.scatter(n_spots, n_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    xlim = pyplot.xlim()
    ylim = pyplot.ylim()
    pyplot.plot([0, max(n_spots)], [0, max(n_spots)], c=red)
    pyplot.xlim(0, xlim[1])
    pyplot.ylim(0, ylim[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('# indexed')
    pyplot.savefig('n_spots_vs_n_indexed.png')
    pyplot.clf()

    pyplot.hexbin(
      n_spots, fraction_indexed, bins='log', cmap=pyplot.cm.jet, gridsize=50)
    pyplot.colorbar()
    #pyplot.scatter(
      #n_spots, fraction_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('Fraction indexed')
    pyplot.savefig('n_spots_vs_fraction_indexed.png')
    pyplot.clf()

    pyplot.hexbin(
      n_indexed, fraction_indexed, bins='log', cmap=pyplot.cm.jet, gridsize=50)
    pyplot.colorbar()
    #pyplot.scatter(
      #n_indexed, fraction_indexed, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# indexed')
    pyplot.ylabel('Fraction indexed')
    pyplot.savefig('n_indexed_vs_fraction_indexed.png')
    pyplot.clf()

    pyplot.hexbin(
      n_spots, n_lattices, bins='log', cmap=pyplot.cm.jet, gridsize=50)
    pyplot.colorbar()
    #pyplot.scatter(
      #n_spots, n_lattices, marker=marker, alpha=alpha, c=blue, lw=lw)
    pyplot.xlim(0, pyplot.xlim()[1])
    pyplot.ylim(0, pyplot.ylim()[1])
    pyplot.xlabel('# spots')
    pyplot.ylabel('# lattices')
    pyplot.savefig('n_spots_vs_n_lattices.png')
    pyplot.clf()

  #pyplot.scatter(
  #  estimated_d_min, d_min_distl_method_1, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.hexbin(estimated_d_min, d_min_distl_method_1, bins='log',
                cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.gca().set_aspect('equal')
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(estimated_d_min), max(d_min_distl_method_1))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('estimated_d_min')
  pyplot.ylabel('d_min_distl_method_1')
  pyplot.savefig('d_min_vs_distl_method_1.png')
  pyplot.clf()

  #pyplot.scatter(
  #  estimated_d_min, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.hexbin(estimated_d_min, d_min_distl_method_2, bins='log',
                cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.gca().set_aspect('equal')
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(estimated_d_min), max(d_min_distl_method_2))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('estimated_d_min')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('d_min_vs_distl_method_2.png')
  pyplot.clf()

  #pyplot.scatter(
  #  d_min_distl_method_1, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.hexbin(d_min_distl_method_1, d_min_distl_method_2, bins='log',
                cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.gca().set_aspect('equal')
  xlim = pyplot.xlim()
  ylim = pyplot.ylim()
  m = max(max(d_min_distl_method_1), max(d_min_distl_method_2))
  pyplot.plot([0, m], [0, m], c=red)
  pyplot.xlim(0, xlim[1])
  pyplot.ylim(0, ylim[1])
  pyplot.xlabel('d_min_distl_method_1')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('distl_method_1_vs_distl_method_2.png')
  pyplot.clf()

  pyplot.hexbin(
    n_spots, estimated_d_min, bins='log', cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.scatter(
    #n_spots, estimated_d_min, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('estimated_d_min')
  pyplot.savefig('n_spots_vs_d_min.png')
  pyplot.clf()

  pyplot.hexbin(
    n_spots, d_min_distl_method_1, bins='log', cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.scatter(
    #n_spots, d_min_distl_method_1, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('d_min_distl_method_1')
  pyplot.savefig('n_spots_vs_distl_method_1.png')
  pyplot.clf()

  pyplot.hexbin(
    n_spots, d_min_distl_method_2, bins='log', cmap=pyplot.cm.jet, gridsize=50)
  pyplot.colorbar()
  #pyplot.scatter(
    #n_spots, d_min_distl_method_2, marker=marker, alpha=alpha, c=blue, lw=lw)
  pyplot.xlim(0, pyplot.xlim()[1])
  pyplot.ylim(0, pyplot.ylim()[1])
  pyplot.xlabel('# spots')
  pyplot.ylabel('d_min_distl_method_2')
  pyplot.savefig('n_spots_vs_distl_method_2.png')
  pyplot.clf()

def spot_count_histogram(n_spots, n_bins=20, filename='n_spots_hist.png', log=False):
  hist = flex.histogram(n_spots.as_double(), n_slots=n_bins)

  blue = '#3498db'

  from matplotlib import pyplot
  pyplot.bar(
    hist.slot_centers().as_numpy_array(),
    hist.slots().as_numpy_array(),
    width=0.75*hist.slot_width(), align='center',
    color=blue, edgecolor=blue, log=log)
  pyplot.savefig(filename)
  pyplot.clf()

def unit_cell_histograms(crystals):
  params = [flex.double() for i in range(6)]
  for cryst in crystals:
    unit_cell = cryst.get_unit_cell().parameters()
    for i in range(6):
      params[i].append(unit_cell[i])

  histograms = []
  for i in range(6):
    histograms.append(flex.histogram(params[i], n_slots=100))

  return histograms

def plot_unit_cell_histograms(crystals):
  histograms = unit_cell_histograms(crystals)

  from matplotlib import pyplot
  f, axes = pyplot.subplots(3, 2, sharey=False, sharex=False)

  blue = '#3498db'
  red = '#e74c3c'

  max_xticks = 5
  xloc = pyplot.MaxNLocator(max_xticks)

  for i, hist in enumerate(histograms):
    col, row = divmod(i, 3)
    ax = axes[row][col]
    ax.bar(
      hist.slot_centers().as_numpy_array(),
      hist.slots().as_numpy_array(),
      width=0.75*hist.slot_width(), align='center',
      color=blue, edgecolor=blue)
    if col == 0:
      ax.set_ylabel('Frequency')
      if row == 2:
        ax.set_xlabel('Unit cell length (Angstrom)')
    elif row == 2:
      ax.set_xlabel('Unit cell angle (degrees)')
  pyplot.savefig('unit_cell_histograms.png')
  pyplot.clf()


if __name__ == '__main__':
  import sys
  run(sys.argv[1:])

