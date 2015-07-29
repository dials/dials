from math import pow
from scitbx import matrix

def determine_miller_ring_sectors(detector, goniometer, s0, hkl_flex, crystal_A):
  crystal_R = matrix.sqr(goniometer.get_fixed_rotation())
  rotation_axis = goniometer.get_rotation_axis()

  from dials.algorithms.spot_prediction import ScanStaticRayPredictor
  from dials.algorithms.spot_prediction import ray_intersection
  from math import radians

  PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI = 1000
  oscillation = (-PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI, PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI)
  rays = ScanStaticRayPredictor(s0, rotation_axis, oscillation)(hkl_flex, crystal_R * crystal_A)
  # ray_intersection could probably be sped up by an is_on_detector() method
  rays = rays.select(ray_intersection(detector, rays))

  divider = radians(10)

  ray_sectors = []
  for s in range(0, 36):
    ray_sectors.append([])
  for (p, m) in zip(rays['phi'], rays['miller_index']):
    ray_sectors[int(p / divider)].append(m)
  return ray_sectors

def geometric_scoring(observations, scorings):
  # Geometric scoring.
  # First observation scores 0.5^0. Second: 0.5^1. Third: 0.5^2. etc.
  # for arbitrary many parallel scoring schemes
  # observations is a list [a']
  # scorings is a dictionary { 'n'->({a':n}, [n]) }
  #    of tuples containing a dictionary mapping a'->int
  #                         and an array of floats representing past observations, onto which the mapping points
  retval = { "scorings": {}, "score": {} }
  total_score = 0.0
  for (name, (scoring, score)) in scorings.iteritems():
    specific_score = 0.0
    score = score[:]
    for o in observations:
      so = scoring[o]
      specific_score += pow(0.5, score[so])
      score[so] += 1
    retval['score'][name] = specific_score
    retval['scorings'][name] = (scoring, score)
    total_score += specific_score
  retval['total'] = total_score
  return retval
