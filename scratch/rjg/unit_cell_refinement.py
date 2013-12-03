from __future__ import division

from dials.array_family import flex
from cctbx import sgtbx, uctbx
import scitbx.lbfgs

def residual(two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  return flex.sum(flex.pow2(two_thetas_obs - two_thetas_calc))

def gradients(
      two_thetas_obs, miller_indices, wavelength, unit_cell, eps=1.e-6):
  result = flex.double()
  for i in xrange(6):
    rs = []
    for signed_eps in [eps, -eps]:
      params_eps = list(unit_cell.parameters())
      params_eps[i] += signed_eps
      rs.append(
        residual(
          two_thetas_obs, miller_indices, wavelength,
          uctbx.unit_cell(params_eps)))
    result.append((rs[0]-rs[1])/(2*eps))
  return result


class refinery(object):

  def __init__(self, two_thetas_obs, miller_indices, wavelength, unit_cell,
               space_group=None):
    self.two_thetas_obs = two_thetas_obs
    self.miller_indices = miller_indices
    #self.unit_cell = unit_cell
    self.space_group = space_group
    self.wavelength = wavelength

    if space_group is not None:
      self.constraints = sgtbx.tensor_rank_2_constraints(
        space_group=self.space_group,reciprocal_space=False)
    else:
      self.constraints = None

    if self.space_group is not None:
      unit_cell = self.space_group.average_unit_cell(unit_cell)

    self.x = flex.double(self.reduce(unit_cell.metrical_matrix()))
    self.minimizer = scitbx.lbfgs.run(target_evaluator=self)

  def reduce(self, metrical_matrix):
    if self.constraints is not None:
      print self.constraints.all_params(
        independent_params=tuple(self.constraints.independent_params(all_params=metrical_matrix)))
      print metrical_matrix
      return self.constraints.independent_params(all_params=metrical_matrix)
    return metrical_matrix

  def enlarge(self,independent):
    # takes reduced number independent parameters, returns 6-parameter metrical matrix
    if self.constraints is not None:
      u_star = self.constraints.all_params(independent_params=tuple(independent))
      assert len(u_star) == 6
      return u_star
    return independent

  def unit_cell(self):
    return uctbx.unit_cell(metrical_matrix=tuple(self.enlarge(self.x)))

  def compute_functional_and_gradients(self):
    unit_cell = self.unit_cell()
    f = residual(
      self.two_thetas_obs, self.miller_indices, self.wavelength, unit_cell)
    g = gradients(
      self.two_thetas_obs, self.miller_indices, self.wavelength, unit_cell)
    g = flex.double(self.reduce(tuple(g)))
    print "functional: %12.6g" % f, "gradient norm: %12.6g" % g.norm()
    return f, g

  def callback_after_step(self, minimizer):
    print "LBFGS step"

def show_fit(two_thetas_obs, miller_indices, wavelength, unit_cell):
  two_thetas_calc = unit_cell.two_theta(miller_indices, wavelength, deg=True)
  for h,o,c in zip(miller_indices, two_thetas_obs, two_thetas_calc):
    print "(%2d, %2d, %2d)" % h, "%6.2f - %6.2f = %6.2f" % (o, c, o-c)
  print


if __name__ == '__main__':
  from cctbx import crystal
  unit_cell_start = uctbx.unit_cell((4.498, 4.498, 7.338, 90, 90, 120))
  crystal_symmetry = crystal.symmetry(
    unit_cell=unit_cell_start,
    space_group_info=sgtbx.space_group_info(number=194))
  ms = crystal_symmetry.build_miller_set(anomalous_flag=False, d_min=1.5)
  ms = ms.sort(by_value="resolution")
  wavelength = 0.97625
  two_thetas = [
    29.3865, 27.2805, 27.1972, 27.1747, 27.364, 27.1991, 27.2159, 27.1775,
    29.4034, 24.9694, 24.9692, 24.9816, 24.9808, 29.4009, 27.1866, 27.1993,
    27.2017, 25.0065, 24.9913, 24.9685, 27.191, 21.0044, 27.1918, 27.1768,
    21.009, 24.9846, 16.2607, 29.4217, 24.982, 27.2008, 15.2354, 15.2498,
    14.3514, 14.352, 21.0061, 16.2576, 24.9884, 21.0089, 15.2646, 16.2589,
    14.346, 14.347, 15.3548, 15.24, 25.0063, 16.2629, 16.266, 15.2613,
    15.2671, 16.2713, 16.2628, 14.3454, 15.2642, 24.9881, 27.4904, 24.975,
    15.2658, 27.2176, 14.3565, 15.2576, 15.2653, 15.2673, 14.3385, 14.355,
    27.2235, 25.0048, 25.0138, 27.1408, 25.0315, 14.3464, 27.2386, 21.0258,
    25.004, 14.3446, 15.2299, 15.2723, 14.3643, 14.3474, 14.3584, 15.2848,
    21.0256, 21.0246, 15.261, 25.0207, 27.2373, 16.2848, 16.2854, 14.3575,
    14.3636, 29.4477, 27.2583, 14.3619, 21.0374, 21.0399, 16.2755, 14.3487,
    14.3618, 14.3608, 15.2829, 27.2497, 15.2715, 15.2699, 16.2646, 16.2786,
    16.2821, 16.2696, 21.0368, 21.0307, 25.0431, 21.0456, 21.0224, 27.2257,
    27.2486, 25.0266, 25.0252, 29.4661, 25.0415, 25.0266, 25.046, 29.4752,
    27.2545, 29.4521, 37.3152, 29.4306, 29.4684, 37.3646, 28.9946, 28.9884,
    29.4736, 29.4737, 30.0142]

  miller_indices = flex.miller_index()
  two_thetas_obs = flex.double()
  for i, two_theta in enumerate(two_thetas):
    d_spacing = uctbx.two_theta_as_d(two_theta, wavelength, deg=True)
    for j, d in enumerate(ms.d_spacings().data()):
      if abs(d - d_spacing) < 0.1:
        miller_indices.append(ms.indices()[j])
        two_thetas_obs.append(two_theta)

  show_fit(
    two_thetas_obs, miller_indices, wavelength, unit_cell_start)

  refined = refinery(
    two_thetas_obs, miller_indices, wavelength, unit_cell_start,
      space_group=crystal_symmetry.space_group())
  print

  show_fit(
    two_thetas_obs, miller_indices, wavelength, refined.unit_cell())

  print "Starting unit cell:", unit_cell_start
  print "Refine unit cell:", refined.unit_cell()
  print
