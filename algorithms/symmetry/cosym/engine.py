"""LBFGS refinement engine for cosym analysis."""
from __future__ import absolute_import, division, print_function

import logging
logger = logging.getLogger(__name__)

import scitbx.lbfgs


class lbfgs_with_curvs(object):
  """Minimise a target function using the LBFGS minimiser.

  Implementation of an LBFGS minimiser using curvature information, according
  to the interface defined by :mod:`scitbx.lbfgs`.

  """

  def __init__(self, target, coords,
               use_curvatures=True,
               termination_params=None):
    """Initialise an lbfgs_with_curvs object.

    Args:
      target (dials.algorithms.target.Target): The target function to minimise.
      coords (scitbx.array_family.flex.double): The starting coordinates for
        minimisation.
      use_curvatures (bool): Whether or not to use curvature information in the
        minimisation. Defaults to True.
      termination_params (scitbx.lbfgs.termination_parameters):
        Override the default termination parameters for the minimisation.

    """
    self.target = target

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
    """Compute the functional, gradients and diagonal.

    Returns:
      tuple: A tuple of the functional, gradients and diagonal, where the
      diagonal is the reciprocal of the curvatures.

    """
    f, g, curvs = self.compute_functional_gradients_and_curvatures()

    # Curvatures of zero will cause a crash, because their inverse is taken.
    assert curvs.all_gt(0.0)
    diags = 1. / curvs
    return f, g, diags

  def curvatures(self):
    """Return the curvatures."""
    return self.target.curvatures(self.x)

  def compute_functional_gradients_and_curvatures(self):
    """Compute the functional, gradients and curvatures.

    Returns:
      tuple: A tuple of the functional, gradients and curvatures.

    """
    self.f, self.g = self.target.compute_functional_and_gradients(self.x)
    self.c = self.curvatures()
    return self.f, self.g, self.c

  def compute_functional_and_gradients(self):
    """Compute the functional and gradients.

    Returns:
      tuple: A tuple of the functional and gradients.

    """
    self.f, self.g = self.target.compute_functional_and_gradients(self.x)
    return self.f, self.g

  def callback_after_step(self, minimizer):
    """Log progress after each successful step of the minimisation."""
    logger.debug('minimization step: f, iter, nfun:')
    logger.debug('%s %s %s' %(self.f, minimizer.iter(), minimizer.nfun()))
