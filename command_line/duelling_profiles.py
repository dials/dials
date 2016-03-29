from __future__ import division

# LIBTBX_SET_DISPATCHER_NAME dev.dials.duelling_profiles

import iotbx.phil

phil_scope = iotbx.phil.parse("""\
  method = *example rt0 flat
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
  rs_node_size = 0.0
    .type = float
  min_isum = None
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

def predict_angles(p0_star, experiment, s0=None):
  '''Predict Ewald sphere crossing angle for RLP x'''

  from scitbx import matrix
  import math

  # Kabsch frame
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
  p_star_m1 = math.sqrt(rho ** 2 - p_star_m3 ** 2)

  p0_star_m1 = p0_star.dot(m1)
  p0_star_m2 = p0_star.dot(m2)
  p0_star_m3 = p0_star.dot(m3)

  cp1 = + p_star_m1 * p0_star_m1 +  p_star_m3 * p0_star_m3
  cp2 = - p_star_m1 * p0_star_m1 +  p_star_m3 * p0_star_m3
  sp1 = + p_star_m1 * p0_star_m3 -  p_star_m3 * p0_star_m1
  sp2 = - p_star_m1 * p0_star_m3 -  p_star_m3 * p0_star_m1

  # FIXME determine order of these (entering, exiting) so usage may be
  # more efficient...

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

  if params.debug:
    print 'angles = %f %f' % angles

  angle = angles[0] if (abs(angles[0] - xyz_mm[2]) <
                        abs(angles[1] - xyz_mm[2])) else angles[1]

  i0 = reflection['intensity.sum.value']# / reflection['dqe']

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
        print '%4d' % data[(j, i)],
      print

  sigma_m = experiment.profile.sigma_m() * d2r
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
    b = random_vector_cone(s0, sigma_b)
    p0 = random_vector_cone(Amat * hkl, sigma_m)
    if params.rs_node_size > 0:
      ns = params.rs_node_size
      import random
      dp0 = matrix.col((random.gauss(0, ns),
                        random.gauss(0, ns),
                        random.gauss(0, ns)))
      p0 += dp0
    angles = predict_angles(p0, experiment, b)
    r = angles[0] if (abs_angle(angles[0], r0) <
                      abs_angle(angles[1], r0)) else angles[1]
    p = p0.rotate(a, r)
    s1 = p + b
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
        print '%4d' % int(patch[(j, i)]),
      print

  print 'Correlation coefficient: %.3f isum: %.1f ' % (
      profile_correlation(data, patch), i0)

  return

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

  nref1 = len(reflections)

  print 'Removed %d invalid reflections, %d remain' % (nref0 - nref1, nref1)

  for j, reflection in enumerate(reflections):
    if ids:
      if j in ids:
        globals()['model_reflection_%s' % method](reflection, experiment, params)
    elif params.id_start and params.id_end:
      if j < params.id_start or j >= params.id_end:
        continue
      print j
      globals()['model_reflection_%s' % method](reflection, experiment, params)

    else:
      print j
      globals()['model_reflection_%s' % method](reflection, experiment, params)

  return

def run(args):
  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_datablocks
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

  main(reflections[0], experiments[0], params)

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
