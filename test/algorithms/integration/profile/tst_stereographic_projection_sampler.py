from __future__ import absolute_import, division, print_function

class Test(object):

  def __init__(self):
    from dxtbx.model.experiment_list import ExperimentListFactory
    import libtbx.load_env
    from os.path import join
    try:
      dials_regression = libtbx.env.dist_path('dials_regression')
    except KeyError:
      print('FAIL: dials_regression not configured')
      exit(0)

    path = join(dials_regression, "centroid_test_data", "experiments.json")

    self.experiments = ExperimentListFactory.from_json_file(path)

  def run(self):

    from dials.array_family import flex
    from math import acos, sqrt, atan, pi, sin, cos
    from scitbx import  matrix
    from dials.algorithms.spot_prediction import ScanStaticRayPredictor

    X = []
    Y = []
    NX = 10
    NY = 10
    for j in range(NX):
      for i in range(NY):
        theta = -pi + i * 2*pi / float(NY)
        phi   = -pi + j * 2*pi / float(NX)
        r = 1
        x = r * sin(theta)*cos(phi)
        y = r * sin(theta)*sin(phi)
        z = -1 + r * cos(theta)
        print(x)
        X.append(x / (1 - z))
        Y.append(y / (1 - z))

    from matplotlib import pylab
    pylab.scatter(X, Y)
    pylab.show()

    # s0 = matrix.col(self.experiments[0].beam.get_s0())
    # m2 = matrix.col(self.experiments[0].goniometer.get_rotation_axis())
    # A  = self.experiments[0].crystal.get_A()
    # dphi = self.experiments[0].scan.get_oscillation_range()
    # s1_list = []

    # zaxis = matrix.col((0, 0, -1))
    # angle = s0.normalize().angle(zaxis, deg=False)
    # axis = s0.normalize().cross(zaxis).normalize()
    # dphi = 0, 2*pi
    # predict = ScanStaticRayPredictor(s0, m2, dphi)

    # for h in range(-20, 20):
    #   for k in range(-20, 20):
    #     for l in range(-20, 20):
    #       if (h,k,l) == (0, 0, 0):
    #         continue
    #       ray = predict((h,k,l), A)
    #       for r in ray:
    #         s1_list.append(r.s1)

    # x = []
    # y = []
    # for i in range(len(s1_list)):
    #   s1 = matrix.col(s1_list[i]).normalize()
    #   s1 = s1.rotate_around_origin(axis, angle, deg=False)
    #   c1 = s1[0] *((1/ (1 - s1[2])))
    #   c2 = s1[1] *((1/ (1 - s1[2])))
    #   # c1 = acos(s1[2])
    #   # c2 = atan(s1[1] / s1[0])
    #   # c1 = s1[0] *sqrt(2/ (1 - s1[2]))
    #   # c2 = s1[1] *sqrt(2/ (1 - s1[2]))
    #   x.append(c1)
    #   y.append(c2)

    # from matplotlib import pylab
    # pylab.scatter(x, y)
    # pylab.show()
    # reflections = flex.reflection_table.from_predictions_multi(self.experiments)

    # print len(reflections)
    # xyz_list = reflections['xyzcal.px']
    # s1_list = reflections['s1']
    # hkl_list = reflections['miller_index']
    # s0 = matrix.col(self.experiments[0].beam.get_s0())
    # m2 = matrix.col(self.experiments[0].goniometer.get_rotation_axis())

    # zaxis = matrix.col((0, 0, -1))
    # angle = s0.normalize().angle(zaxis, deg=False)
    # axis = s0.normalize().cross(zaxis).normalize()

    # phi_list = []
    # theta_list = []
    # dx = []
    # dy = []
    # dz = []
    # x, y, z = zip(*list(xyz_list))

    # for i, s1 in enumerate(s1_list):
    #   s1 = matrix.col(s1).normalize()
    #   s1 = s1.rotate_around_origin(axis, angle, deg=False)
    #   c3 = xyz_list[i][2]
    #   c1 = s1[0] *sqrt(2/ (1 - s1[2]))
    #   c2 = s1[1] *sqrt(2/ (1 - s1[2]))
    #   dx.append(c1)
    #   dy.append(c2)
    #   dz.append(c3)


    #   # theta = acos(s1[2] / s1.length())
    #   # phi = atan(s1[1] / s1[0])
    #   # theta_list.append(theta)
    #   # phi_list.append(phi)
    #   # dx.append(xyz_list[i][0])
    #   # dy.append(xyz_list[i][1])

    # # from dials.algorithms.spatial_indexing import make_spatial_index

    # # coords = flex.vec3_double([(x,y,z) for x,y,z in zip(dx,dy,dz)])

    # # index = make_spatial_index(coords, 50)
    # # print index
    # # print dir(index)

    # from matplotlib import pylab
    # # pylab.scatter(theta_list, phi_list)
    # # pylab.show()
    # # pylab.scatter(dx, dy)
    # # pylab.scatter(x, y)
    # pylab.scatter(dx, dy)
    # # pylab.hexbin(theta_list, phi_list, C=dx, gridsize=30)
    # pylab.show()

if __name__ == '__main__':

  test = Test()
  test.run()
