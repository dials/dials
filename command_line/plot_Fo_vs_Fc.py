#!/usr/bin/env python

# LIBTBX_SET_DISPATCHER_NAME dev.dials.plot_Fo_vs_Fc

"""
Create a plot of Fo vs Fc similar to that shown by Figure 6 in
https://doi.org/10.1107/S2059798317010348

Usage: dev.dials.plot_Fo_vs_Fc hklin=refined.mtz
"""

from __future__ import division, print_function, absolute_import
import sys
from dials.util import Sorry
from dials.util.options import OptionParser
#from libtbx.table_utils import simple_table
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator
from iotbx import mtz
from scitbx.array_family import flex
from scitbx.lstbx import normal_eqns, normal_eqns_solving
from math import sqrt

class HyperbolaFit(normal_eqns.non_linear_ls, normal_eqns.non_linear_ls_mixin):
  """Fit the function y = sqrt(x^2 + a^2) by non-linear regression. There is
  just one parameter, a^2."""

  # Initial guess for the value of a^2
  a_sq0 = flex.double([1000])

  def __init__(self, x, y):
    super(HyperbolaFit, self).__init__(n_parameters=1)
    self.x = x
    self.y = y
    self.n_data = len(self.x)
    assert len(self.y) == self.n_data
    self.restart()

  def restart(self):
    self.param = self.a_sq0.deep_copy()
    self.old_param = None

  def parameter_vector_norm(self):
    return self.param.norm()

  def build_up(self, objective_only=False):
    a_sq = self.param[0]
    model_y = flex.sqrt(flex.pow2(self.x) + a_sq)
    residuals = model_y - self.y

    self.reset()
    if objective_only:
      self.add_residuals(residuals, weights=None)
    else:
      dy_dp = 0.5/model_y
      jacobian = flex.double(flex.grid(self.n_data, 1))
      jacobian.matrix_paste_column_in_place(dy_dp, 0)
      self.add_equations(residuals, jacobian, weights=None)

  def step_forward(self):
    self.old_param = self.param.deep_copy()
    self.param += self.step()

  def step_backward(self):
    assert self.old_param is not None
    self.param, self.old_param = self.old_param, None

  def goodness_of_fit(self):
    """Calculate various goodness of fit metrics (assumes fit has been
    performed already)"""
    a_sq = self.param[0]
    model_y = flex.sqrt(flex.pow2(self.x) + a_sq)
    resid = model_y - self.y
    resid2 = flex.pow2(resid)

    sse = flex.sum(resid2)
    sst = flex.sum(flex.pow2(model_y - flex.mean(model_y)))
    r_sq = 1 - sse/sst
    rmse = sqrt(sse / (self.n_data - 1))

    return {"SSE": sse, "R-square":r_sq, "RMSE":rmse}

class Script(object):
  '''A class for running the script.'''

  def __init__(self):
    '''Initialise the script.'''
    from libtbx.phil import parse
    import libtbx.load_env

    # The phil scope
    phil_scope = parse('''
      hklin = None
        .type = path
        .help = "MTZ file (containing observed and calculated structure factors)"

      Fo = F
        .type = str
        .help = "MTZ column name for Fobs"

      Fc = FC_ALL
        .type = str
        .help = "MTZ column name for Fcalc (FC_ALL from Refmac includes the"
                "bulk solvent contribution)"

      max_Fc = 300
        .type = float
        .help = "Set plot limits to display data up to this value of Fc"

      plot_filename = Fo_vs_Fc.pdf
        .type = str
        .help = "Filename for plot"

      fit_hyperbola = True
        .type = bool
        .help = "Calculate and show the fit of a hyperbolic function given by"
                "|Fo|^2 = |Fc|^2 + |Fe|^2, where |Fe| describes the error term"
                "containing information about dynamic scattering and other"
                "effects"

      show_y_eq_x = True
        .type = bool
        .help = "Plot y=x as a dashed line"
    ''', process_includes=True)

    # The script usage
    usage = ("usage: dev.dials.plot_Fo_vs_Fc hklin=refined.mtz")

    # Create the parser
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      epilog=__doc__)

    self.model_fit=None

    return

  def _extract_data_from_mtz(self):
    try:
      m = mtz.object(self.params.hklin)
    except RuntimeError:
      raise Sorry('Could not read {0}'.format(self.params.hklin))

    mad = m.as_miller_arrays_dict()
    mad = {k[-1]: v for (k,v) in mad.items()}
    fobs = mad.get(self.params.Fo)
    fc = mad.get(self.params.Fc)

    if [fobs, fc].count(None) > 0:
      raise Sorry('Columns {0} not found in available labels: {1}'.format(
        ", ".join([self.params.Fo, self.params.Fc]),
        ", ".join(m.column_labels())))

    # Find common reflections (some fobs might be missing)
    fobs, fc = fobs.common_sets(fc)

    self.fobs = fobs.data()
    self.fc = fc.amplitudes().data()

    return

  def _plot(self):
    fig = plt.figure()
    ax = fig.add_subplot(111)
    minor_loc = MultipleLocator(10)
    ax.yaxis.set_minor_locator(minor_loc)
    ax.xaxis.set_minor_locator(minor_loc)
    ax.grid(True, which='minor')
    ax.set_axisbelow(True)
    ax.set_aspect('equal')
    ax.set_xlabel(r'$F_c$')
    ax.set_ylabel(r'$F_o$')
    ax.scatter(self.fc, self.fobs, s=1, c="indianred")

    if self.params.max_Fc:
      ax.set_xlim((0, self.params.max_Fc))
      ax.set_ylim((0, self.params.max_Fc))

    if self.params.show_y_eq_x:
      ax.plot(ax.get_xlim(), ax.get_ylim(), ls="--", c="0.0", linewidth=0.8)

    if self.model_fit:
      x = flex.double_range(0, int(ax.get_xlim()[1]))
      y = self.model_fit(x)
      ax.plot(x, y, c="0.0", linewidth=0.8)

    print("Saving plot to {0}".format(self.params.plot_filename))
    plt.savefig(self.params.plot_filename)

  def run(self):
    '''Execute the script.'''

    # Parse the command line
    self.params, options = self.parser.parse_args(show_diff_phil=True)

    if self.params.hklin is None:
      self.parser.print_help()
      sys.exit()

    self._extract_data_from_mtz()

    if self.params.fit_hyperbola:
      # fit by NLLS Levenberg Marquardt algorithm
      hyperbola_fit = HyperbolaFit(self.fc, self.fobs)
      hyperbola_fit.restart()
      iterations = normal_eqns_solving.levenberg_marquardt_iterations(
        hyperbola_fit,
        track_all=True,
        gradient_threshold=1e-8,
        step_threshold=1e-8,
        tau=1e-4,
        n_max_iterations=200)
      intercept = hyperbola_fit.param[0]

      print("Model fit described by the formula: |Fo|^2 = sqrt(|Fc|^2 + |Fe|^2)")
      print("where |Fe| = {:.5f}\n".format(sqrt(intercept)))

      print("Goodness of fit:")
      gof = hyperbola_fit.goodness_of_fit()
      print("SSE: {:.5g}".format(gof['SSE']))
      print("R-square: {:.5f}".format(gof['R-square']))
      print("RMSE: {:.2f}".format(gof['RMSE']))
      print()

      # Set the model_fit function using the determined intercept
      def hyperbola(x, c):
        return flex.sqrt(flex.pow2(x) + c)
      from functools import partial
      self.model_fit = partial(hyperbola, c=intercept)

    if self.params.plot_filename:
      self._plot()

    return

if __name__ == '__main__':
  from dials.util import halraiser

  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
