from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

import scitbx.lbfgs


class lbfgs_with_curvs(object):
  def __init__(self, target, coords,
               animate=False,
               save_intermediate_plots=False,
               use_curvatures=True,
               termination_params=None):
    self.target = target
    self.animate = animate
    self.save_intermediate_plots = save_intermediate_plots

    self.dim = len(coords)
    self.x = coords

    if use_curvatures:
      self.diag_mode = "always"
    else:
      self.diag_mode = None

    self.minimizer = scitbx.lbfgs.run(
      target_evaluator=self,
      termination_params=termination_params
    )

  def compute_functional_gradients_diag(self):
    f, g, curvs = self.compute_functional_gradients_and_curvatures()

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert curvs.all_gt(0.0)
    diags = 1. / curvs
    return f, g, diags

  def curvatures(self):
    return self.target.curvatures(self.x)

  def compute_functional_gradients_and_curvatures(self):
    self.f, self.g = self.target.compute_functional_and_gradients(self.x)
    self.c = self.curvatures()
    return self.f, self.g, self.c

  def compute_functional_and_gradients(self):
    self.f, self.g = self.target.compute_functional_and_gradients(self.x)
    return self.f, self.g

  def callback_after_step(self, minimizer):
    logger.debug('minimization step: f, iter, nfun:')
    logger.debug('%s %s %s' %(self.f, minimizer.iter(), minimizer.nfun()))

    if self.animate:
      from matplotlib import pyplot as plt
      NN = self.x.size() // self.dim
      coord_x = self.x[0:NN]
      coord_y = self.x[NN:2*NN]
      plt.clf()
      plt.ion()
      plt.scatter(coord_x, coord_y, c="r", marker='+', s=3)
      plt.axes().set_aspect("equal")
      plt.xlim(-1, 1)
      plt.ylim(-1, 1)
      plt.pause(0.005)

    if self.save_intermediate_plots:
      from dials.algorithms.symmetry import cosym
      NN = self.x.size() // self.dim
      coord_x = self.x[0:NN]
      coord_y = self.x[NN:2*NN]
      cosym.plot((coord_x, coord_y), plot_name='xy_step_%02i.png' %minimizer.iter())

