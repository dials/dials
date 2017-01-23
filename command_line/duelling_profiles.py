from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.duelling_profiles

import iotbx.phil

phil_scope = iotbx.phil.parse("""\
  method = *example rt0 flat predict
    .type = choice
  id = None
    .type = int
    .multiple = True
  id_start = None
    .type = int
  id_end = None
    .type = int
  scale = 1.0
    .type = float
  show = False
    .type = bool
  physics = False
    .type = bool
  rs_node_size = 0.0
    .type = float
  min_isum = None
    .type = float
  num = -1
    .type = int
  seed = -1
    .type = int
  sigma_m = 0
    .type = float
  sigma_b = 0
    .type = float
  sigma_l = 0
    .type = float
  sigma_cell = 0
    .type = float
  debug = False
    .type = bool
""", process_includes=True)

help_message = '''

Examples::

  dev.dials.duelling_profiles experiments.json integrated.pickle

'''

def model_background(shoebox, mean_bg):
  from scitbx.random import variate, poisson_distribution
  dz, dy, dx = shoebox.focus()
  g = variate(poisson_distribution(mean = mean_bg))
  for k in range(dz):
    for j in range(dy):
      for i in range(dx):
        shoebox[k, j, i] += g.next()
  return

def random_vector_2D(vector, sd):
  '''Form random orthonormal basis including vector, rotate vector by random
  amount sd in both normal directions.'''
  import random
  o0 = vector.ortho()
  o1 = vector.cross(o0)
  return vector.rotate(o0, random.gauss(0, sd)).rotate(o1, random.gauss(0, sd))

def random_vector_cone(vector, sd):
  '''Pick orthogonal axis to vector, rotate by random angle about vector,
  rotate vector about this by sd in radians.'''
  import random
  import math
  o0 = vector.ortho()
  o1 = vector.rotate(o0, random.random() * 2.0 * math.pi)
  return vector.rotate(o1, random.gauss(0, sd))

def model_reflection_example(reflection, experiment, params):
  hkl = reflection['miller_index']
  i0 = reflection['intensity.sum.value'] / reflection['dqe']
  s1 = reflection['s1']
  xyz = reflection['xyzcal.px']
  pixels = reflection['shoebox']
  mean_bg = reflection['background.mean']
  crystal = experiment.crystal
  profile = experiment.profile
  Amat = crystal.get_A_at_scan_point(int(xyz[2]))
  return

def trace_path(x0, y0, v, f, s, t0, p):
  '''x0, y0 impact position in mm, v, f, s are vectors of the ray and fast
  and slow directions in lab frame, t0 is thickness in mm and p pixel size in
  mm'''

  import math

  n = f.cross(s)
  a = v.angle(n)
  t = t0 / math.cos(a)

  dx = t0 * v.dot(f) / v.dot(n)
  dy = t0 * v.dot(s) / v.dot(n)
  x1 = x0 + dx
  y1 = y0 + dy

  nx0 = x0 / p
  nx1 = x1 / p
  ny0 = y0 / p
  ny1 = y1 / p

  pixels = []

  if abs(dy) > abs(dx):
    if dy > 0:
      start, end, step = int(ny0), int(ny1) + 1, 1
      d0, d1 = 0, 1
    else:
      start, end, step = int(ny0), int(ny1) - 1, -1
      d0, d1 = -1, 0

    m = dx / dy
    s = int(round(m / abs(m))) * step
    pixels = []
    for j in range(start, end, step):
      if j == start:
        l0 = ny0
      else:
        l0 = j + d0 * step
      if j == (end - step):
        l1 = ny1
      else:
        l1 = j + d1 * step

      m0 = nx0 + (l0 - ny0) * m
      m1 = nx0 + (l1 - ny0) * m
      c = int(m0) if s < 0 else int(m1)
      if int(m1) != int(m0):
        l2 = l0 + (c - m0) / m
        _s = p * math.sqrt((l2 - l0) ** 2 + (c - m0) ** 2) / math.sin(a)
        pixels.append((int(m0), int(l0), _s))
        _s = p * math.sqrt((l1 - l2) ** 2 + (m1 - c) ** 2) / math.sin(a)
        pixels.append((int(m1), int(l0), _s))
      else:
        _s = p * math.sqrt((l1 - l0) ** 2 + (m1 - m0) ** 2) / math.sin(a)
        pixels.append((int(c), int(l0), _s))

  else:
    if dx > 0:
      start, end, step = int(nx0), int(nx1) + 1, 1
      d0, d1 = 0, 1
    else:
      start, end, step = int(nx0), int(nx1) - 1, -1
      d0, d1 = -1, 0
    m = dy / dx
    s = int(round(m / abs(m))) * step
    for j in range(start, end, step):
      if j == start:
        l0 = nx0
      else:
        l0 = j + d0 * step
      if j == (end - step):
        l1 = nx1
      else:
        l1 = j + d1 * step

      m0 = ny0 + (l0 - nx0) * m
      m1 = ny0 + (l1 - nx0) * m
      c = int(m0) if s < 0 else int(m1)
      if int(m1) != int(m0):
        l2 = l0 + (c - m0) / m
        _s = p * math.sqrt((l2 - l0) ** 2 + (c - m0) ** 2) / math.sin(a)
        pixels.append((int(l0), int(m0), _s))
        _s = p * math.sqrt((l1 - l2) ** 2 + (m1 - c) ** 2) / math.sin(a)
        pixels.append((int(l0), int(m1), _s))
      else:
        _s = p * math.sqrt((l1 - l0) ** 2 + (m1 - m0) ** 2) / math.sin(a)
        pixels.append((int(l0), int(c), _s))
  return pixels


def model_path_through_sensor(detector, reflection, s1, patch, scale):
  '''Model the passage of the ray s1 through the detector, depositing
  fractional counts in patch as we go.'''

  import math
  from scitbx import matrix

  p, xy = detector.get_ray_intersection(s1)
  x0, y0 = xy

  panel = detector[p]

  mu = panel.get_mu()
  t0 = panel.get_thickness()
  pixel = panel.get_pixel_size()
  f = matrix.col(panel.get_fast_axis())
  s = matrix.col(panel.get_slow_axis())
  n = f.cross(s)
  t = t0 / math.cos(s1.angle(n))

  v = s1.normalize()

  pixels = trace_path(x0, y0, v, f, s, t0, pixel[0])

  photon = scale

  bbox = reflection['bbox']

  for x, y, l in pixels:
    deposit = photon * (1 - math.exp(-mu * l))
    photon -= deposit
    if x < bbox[0] or x >= bbox[1]:
      continue
    if y < bbox[2] or y >= bbox[3]:
      continue
    x -= bbox[0]
    y -= bbox[2]
    patch[(y, x)] += deposit

  return scale - photon

def model_reflection_predict(reflection, experiment, params):
  import math
  from scitbx import matrix
  from dials.array_family import flex

  d2r = math.pi / 180.0

  hkl = reflection['miller_index']
  xyz = reflection['xyzcal.px']
  xyz_mm = reflection['xyzcal.mm']

  if params.debug:
    print 'hkl = %d %d %d' % hkl
    print 'xyz px = %f %f %f' % xyz
    print 'xyz mm = %f %f %f' % xyz_mm

  Amat = matrix.sqr(experiment.crystal.get_A_at_scan_point(int(xyz[2])))
  p0_star = Amat * hkl

  angles = predict_angles(p0_star, experiment)

  assert(angles)

  if params.debug:
    print 'angles = %f %f' % angles

  angle = angles[0] if reflection['entering'] else angles[1]

  if abs_angle(angle, xyz_mm[2]) > 1.0e-3:
    raise RuntimeError, '%f %f' % (angle, xyz_mm[2])

  return

def predict_angles(p0_star, experiment, s0=None):
  '''Predict Ewald sphere crossing angle for RLP x, returned in order
  (entering, exiting).'''

  from scitbx import matrix
  import math

  a = matrix.col(experiment.goniometer.get_rotation_axis())
  if s0 is None:
    b = matrix.col(experiment.beam.get_s0())
  else:
    b = s0

  m2 = a.normalize()
  m1 = m2.cross(b.normalize())
  m3 = m1.cross(m2)

  p0_sqr = p0_star.dot(p0_star)
  rho = math.sqrt(p0_sqr - p0_star.dot(m2) ** 2)
  p_star_m3 = (-0.5 * p0_sqr - p0_star.dot(m2) * b.dot(m2)) / b.dot(m3)
  p_star_m2 = p0_star.dot(m2)
  if rho ** 2 < p_star_m3 ** 2:
    return None
  p_star_m1 = math.sqrt(rho ** 2 - p_star_m3 ** 2)

  p0_star_m1 = p0_star.dot(m1)
  p0_star_m2 = p0_star.dot(m2)
  p0_star_m3 = p0_star.dot(m3)

  cp1 = + p_star_m1 * p0_star_m1 +  p_star_m3 * p0_star_m3
  cp2 = - p_star_m1 * p0_star_m1 +  p_star_m3 * p0_star_m3
  sp1 = + p_star_m1 * p0_star_m3 -  p_star_m3 * p0_star_m1
  sp2 = - p_star_m1 * p0_star_m3 -  p_star_m3 * p0_star_m1

  return math.atan2(sp1, cp1), math.atan2(sp2, cp2)

def profile_correlation(data, model):
  '''Compute CC between reflection profiles data and model.'''

  from dials.array_family import flex

  assert(data.focus() == model.focus())

  sel = data.as_1d() > 0

  model = model.as_1d().select(sel)
  data = data.as_1d().select(sel).as_double()

  assert(len(data) == len(model))

  correlation = flex.linear_correlation(data, model)
  return correlation.coefficient()

def abs_angle(a, b):
  import math
  return abs(math.atan2(math.sin(b - a), math.cos(b - a)))

def model_reflection_rt0(reflection, experiment, params):
  import math
  import random
  from scitbx import matrix
  from dials.array_family import flex

  d2r = math.pi / 180.0

  hkl = reflection['miller_index']
  xyz = reflection['xyzcal.px']
  xyz_mm = reflection['xyzcal.mm']

  if params.debug:
    print 'hkl = %d %d %d' % hkl
    print 'xyz px = %f %f %f' % xyz
    print 'xyz mm = %f %f %f' % xyz_mm
    if reflection['entering']:
      print 'entering'
    else:
      print 'exiting'

  Amat = matrix.sqr(experiment.crystal.get_A_at_scan_point(int(xyz[2])))
  p0_star = Amat * hkl

  angles = predict_angles(p0_star, experiment)

  assert(angles)

  if params.debug:
    print 'angles = %f %f' % angles

  angle = angles[0] if (abs(angles[0] - xyz_mm[2]) <
                        abs(angles[1] - xyz_mm[2])) else angles[1]

  p = experiment.detector[reflection['panel']]
  n = matrix.col(p.get_normal())
  s1 = matrix.col(reflection['s1'])
  t = p.get_thickness() / math.cos(s1.angle(n))
  if params.debug:
    print 'dqe = %f' % reflection['dqe']

  if params.physics:
    i0 = reflection['intensity.sum.value'] / reflection['dqe']
  else:
    i0 = reflection['intensity.sum.value']

  if params.min_isum:
    if i0 < params.min_isum:
      return

  s1 = reflection['s1']
  a = matrix.col(experiment.goniometer.get_rotation_axis())
  s0 = matrix.col(experiment.beam.get_s0())

  if params.debug:
    print 's1 = %f %f %f' % s1

  pixels = reflection['shoebox']
  pixels.flatten()
  data = pixels.data
  dz, dy, dx = data.focus()

  # since now 2D data
  data.reshape(flex.grid(dy, dx))

  if params.show:
    print 'Observed reflection (flattened in Z):'
    print
    for j in range(dy):
      for i in range(dx):
        print '%5d' % data[(j, i)],
      print

  if params.sigma_m > 0:
    sigma_m = params.sigma_m * d2r
  else:
    sigma_m = experiment.profile.sigma_m() * d2r

  if params.sigma_b > 0:
    sigma_b = params.sigma_b * d2r
  else:
    sigma_b = experiment.profile.sigma_b() * d2r

  r0 = xyz_mm[2]

  detector = experiment.detector

  patch = flex.double(dy * dx, 0)
  patch.reshape(flex.grid(dy, dx))

  bbox = reflection['bbox']

  scale = params.scale
  if params.show:
    print '%d rays' % (int(round(i0 * scale)))
  for i in range(int(round(i0 * scale))):
    if params.sigma_l:
      l_scale = random.gauss(1, params.sigma_l)
      b = random_vector_cone(s0 * l_scale, sigma_b)
    else:
      b = random_vector_cone(s0, sigma_b)
    if params.sigma_cell:
      cell_scale = random.gauss(1, params.sigma_cell)
      p0 = random_vector_cone(cell_scale * Amat * hkl, sigma_m)
    else:
      p0 = random_vector_cone(Amat * hkl, sigma_m)
    if params.rs_node_size > 0:
      ns = params.rs_node_size
      import random
      dp0 = matrix.col((random.gauss(0, ns),
                        random.gauss(0, ns),
                        random.gauss(0, ns)))
      p0 += dp0
    angles = predict_angles(p0, experiment, b)
    if angles is None:
      # scattered ray ended up in blind region
      continue
    r = angles[0] if reflection['entering'] else angles[1]
    p = p0.rotate(a, r)
    s1 = p + b

    if params.physics:
      model_path_through_sensor(detector, reflection, s1, patch, scale)

    else:
      panel, xy = detector.get_ray_intersection(s1)

      # FIXME DO NOT USE THIS FUNCTION EVENTUALLY...
      x, y = detector[panel].millimeter_to_pixel(xy)
      if x < bbox[0] or x >= bbox[1]:
        continue
      if y < bbox[2] or y >= bbox[3]:
        continue
      x -= bbox[0]
      y -= bbox[2]
      # FIXME in here try to work out probability distribution along path
      # length through the detector sensitive surface i.e. deposit fractional
      # counts along pixels (and allow for DQE i.e. photon passing right through
      # the detector)
      patch[(int(y), int(x))] += 1.0 / scale

  if params.show:
    print 'Simulated reflection (flattened in Z):'
    print
    for j in range(dy):
      for i in range(dx):
        print '%5d' % int(patch[(j, i)]),
      print


  cc = profile_correlation(data, patch)
  print 'Correlation coefficient: %.3f isum: %.1f ' % (cc, i0)

  return cc

def model_reflection_flat(reflection, experiment, params):
  pixels = reflection['shoebox']
  pixels.flatten()
  data = pixels.data
  dz, dy, dx = data.focus()
  return

def main(reflections, experiment, params):
  nref0 = len(reflections)

  method = params.method
  ids = params.id

  if 'intensity.prf.variance' in reflections:
    selection = reflections.get_flags(
      reflections.flags.integrated,
      all=True)
  else:
    selection = reflections.get_flags(
      reflections.flags.integrated_sum)
  reflections = reflections.select(selection)

  # filter according to rules

  if params.min_isum:
    selection = reflections['intensity.sum.value'] > params.min_isum
    reflections = reflections.select(selection)

  if params.num > len(reflections):
    raise RuntimeError, 'you asked for too many reflections sorry'

  if params.seed > 0 and params.num > 0:
    import random
    from dials.array_family import flex
    random.seed(params.seed)
    selected = flex.bool(len(reflections), False)
    while len(selected.iselection()) < params.num:
      selected[random.randint(0, len(reflections))] = True
    reflections = reflections.select(selected)

  nref1 = len(reflections)

  print 'Removed %d reflections, %d remain' % (nref0 - nref1, nref1)

  results = []

  for j, reflection in enumerate(reflections):
    result = None
    if ids:
      if j in ids:
        result = globals()['model_reflection_%s' % method](reflection, experiment, params)
    elif params.id_start and params.id_end:
      if j < params.id_start or j >= params.id_end:
        continue
      print j
      result = globals()['model_reflection_%s' % method](reflection, experiment, params)

    else:
      print j
      result = globals()['model_reflection_%s' % method](reflection, experiment, params)
    if result is not None:
      results.append(result)

  return results

def run(args):
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  import libtbx.load_env

  usage = "%s [options] integrated.pickle experiments.json" % (
    libtbx.env.dispatcher_name)

  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True,
    check_format=False,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=False)
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  if len(experiments) != 1 or len(reflections) != 1:
    parser.print_help()
    exit()

  if not 'shoebox' in reflections[0]:
    print 'Please add shoeboxes to reflection pickle'
    exit()

  results = main(reflections[0], experiments[0], params)

  if results:
    print 'mean result: %.3f' % (sum(results) / len(results))

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
