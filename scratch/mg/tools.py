    def determine_miller_ring_sectors(detector, goniometer, s0, possible_hkl, crystal_A):
      crystal_R = matrix.sqr(goniometer.get_fixed_rotation())
      rotation_axis = goniometer.get_rotation_axis()

      from dials.algorithms.spot_prediction import ScanStaticRayPredictor
      from dials.algorithms.spot_prediction import ray_intersection
      from math import radians

      PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI = 1000
      oscillation = (-PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI, PRACTICALLY_INFINITY_BUT_DEFINITELY_LARGER_THAN_2PI)
      rays = ScanStaticRayPredictor(s0, rotation_axis, oscillation)(flex.miller_index(possible_hkl), crystal_R * crystal_A)
      # ray_intersection could probably be sped up by an is_on_detector() method
      rays = rays.select(ray_intersection(detector, rays))

      divider = radians(10)

      ray_sectors = []
      for s in range(0, 36):
        ray_sectors.append([])
      for (p, m) in zip(rays['phi'], rays['miller_index']):
        ray_sectors[int(p / divider)].append(m)
      return ray_sectors
