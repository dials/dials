from __future__ import division
#
#  Copyright (C) (2016) Lawrence Berkeley National Laboratory
#
#  Author: Nicholas Sauter.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#
import libtbx
from libtbx.utils import Sorry
from scitbx.array_family import flex
from logging import debug
from engine import TARGET_ACHIEVED,RMSD_CONVERGED,STEP_TOO_SMALL
from engine import MAX_ITERATIONS,MAX_TRIAL_ITERATIONS,DOF_TOO_LOW

try:
  from scitbx.examples.bevington import non_linear_ls_eigen_wrapper
except ImportError,e:
  raise Sorry("""Eigen package is not available.  Please untar the Eigen source package
     (http://eigen.tuxfamily.org) and place a link to it (eigen--> Eigen source dir) in
     the modules directory of your developer install; then recompile.
""")

from engine import AdaptLstbx as AdaptLstbxBase

class AdaptLstbxSparse(AdaptLstbxBase,non_linear_ls_eigen_wrapper):
  """Adapt the base class for Eigen"""

  def __init__(self, target, prediction_parameterisation, log=None,
               verbosity = 0, track_step = False, track_gradient = False,
               track_parameter_correlation = False,
               track_out_of_sample_rmsd = False, max_iterations = None):

    AdaptLstbxBase.__init__(self, target, prediction_parameterisation,
             log=log, track_step=track_step,
             track_gradient=track_gradient,
             track_parameter_correlation=track_parameter_correlation,
             track_out_of_sample_rmsd=track_out_of_sample_rmsd,
             max_iterations=max_iterations)

    non_linear_ls_eigen_wrapper.__init__(self, n_parameters = len(self._parameters))

  def set_nproc(self, *args, **kwargs):
    # Multiprocessing is not implemented for sparse matrix LevMar
    # XXX possible future path: pickle support for sparse::matrix, then divide up the Jacobian
    self._nproc = 1
    return

  def finalise(self):
    """in contrast to dense matrix (traditional) LevMar, sparse matrix assumes
       that the matrix is extremely large and not easily inverted.  Therefore,
       no attempt here to calculate esd's based on the variance covariance matrix."""
    self.show_eigen_summary()
    return

from engine import GaussNewtonIterations as GaussNewtonIterationsBase
class GaussNewtonIterations(AdaptLstbxSparse, GaussNewtonIterationsBase):
  """Refinery implementation, using lstbx Gauss Newton iterations"""

  def __init__(self, target, prediction_parameterisation, log=None,
               verbosity=0, track_step=False, track_gradient=False,
               track_parameter_correlation=False,
               track_out_of_sample_rmsd=False,
               max_iterations=20, **kwds):

    AdaptLstbxSparse.__init__(self, target, prediction_parameterisation,
             log=log, verbosity=verbosity, track_step=track_step,
             track_gradient=track_gradient,
             track_parameter_correlation=track_parameter_correlation,
             track_out_of_sample_rmsd=track_out_of_sample_rmsd,
             max_iterations=max_iterations)

    # add an attribute to the journal
    self.history.add_column("reduced_chi_squared")#flex.double()

    # adopt any overrides of the defaults above
    libtbx.adopt_optional_init_args(self, kwds)

from engine import LevenbergMarquardtIterations
class SparseLevenbergMarquardtIterations(GaussNewtonIterations,LevenbergMarquardtIterations):
  """Levenberg Marquardt with Sparse matrix algebra"""

  def run(self):

    # add an attribute to the journal
    self.history.add_column("mu")
    self.history.add_column("nu")

    #FIXME need a much neater way of doing this stuff through
    #inheritance
    # set max iterations if not already.
    if self._max_iterations is None:
      self._max_iterations = 20

    self.n_iterations = 0
    nu = 2
    self.build_up()

    # return early if refinement is not possible
    if self.dof < 1:
      self.history.reason_for_termination = DOF_TOO_LOW
      return

    a_diag = self.get_normal_matrix_diagonal()
    self.mu = self.tau * flex.max(a_diag)

    while True:

      # set functional and gradients for the step
      self._f = self.objective()
      self._g = -self.opposite_of_gradient()

      # cache some items for the journal prior to solve
      pvn = self.parameter_vector_norm()
      gn = self.opposite_of_gradient().norm_inf()

      self.add_constant_to_diagonal(self.mu)

      # solve the normal equations
      self.solve()

      # standard journalling
      self.update_journal()
      debug("Step %d", self.history.get_nrows() - 1)

      # add cached items to the journal
      self.history.set_last_cell("parameter_vector_norm", pvn)
      self.history.set_last_cell("gradient_norm", gn)

      # extra journalling post solve
      self.history.set_last_cell("mu", self.mu)
      self.history.set_last_cell("nu", nu)
      if self.history.has_key("solution"):
        self.history.set_last_cell("solution", self.actual.step().deep_copy())
      self.history.set_last_cell("solution_norm", self.step().norm())
      self.history.set_last_cell("reduced_chi_squared", self.chi_sq())

      # test termination criteria before taking the next forward step
      if self.had_too_small_a_step():
        self.history.reason_for_termination = STEP_TOO_SMALL
        break
      if self.test_for_termination():
        self.history.reason_for_termination = TARGET_ACHIEVED
        break
      if self.test_rmsd_convergence():
        self.history.reason_for_termination = RMSD_CONVERGED
        break
      if self.n_iterations == self._max_iterations:
        self.history.reason_for_termination = MAX_ITERATIONS
        break

      h = self.step()
      expected_decrease = 0.5*h.dot(self.mu*h - self._g)
      self.step_forward()
      self.n_iterations += 1
      self.build_up(objective_only=True)
      objective_new = self.objective()
      print "%5d %18.4f"%(self.n_iterations,objective_new), "%12.7f"%(self.mu)
      actual_decrease = self._f - objective_new
      rho = actual_decrease/expected_decrease
      if rho > 0:
        self.mu *= max(1/3, 1 - (2*rho - 1)**3)
        nu = 2
      else:
        self.step_backward()
        self.history.del_last_row()
        if nu >= 8192:
          self.history.reason_for_termination = MAX_TRIAL_ITERATIONS
          break
        self.mu *= nu
        nu *= 2

      # prepare for next step
      self.build_up()

    self.finalise()

    return
