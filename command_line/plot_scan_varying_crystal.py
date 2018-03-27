#!/usr/bin/env python
#
# plot_scan_varying_crystal.py
#
#  Copyright (C) (2014) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David G. Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

from __future__ import absolute_import, division, print_function

import os
import matplotlib
matplotlib.use('Agg')

from dials.algorithms.refinement.rotation_decomposition import \
  solve_r3_rotation_for_angles_given_axes
from dials.command_line.analyse_output import ensure_directory

from libtbx.utils import Sorry
from libtbx.phil import parse
phil_scope = parse('''
  output {
    directory = .
      .type = str
      .help = "The directory to store the results"

    format = *png pdf
      .type = choice

    debug = False
      .help = "print tables of values that will be plotted"
      .type = bool
      .expert_level = 1
  }

  orientation_decomposition
    .help = "Options determining how the orientation matrix"
            "decomposition is done. The axes about which to decompose"
            "the matrix into three rotations are chosen here, as well"
            "as whether the rotations are relative to the reference"
            "orientation, taken from the static crystal model"
  {
    e1 = 1. 0. 0.
      .type = floats(size = 3)

    e2 = 0. 1. 0.
      .type = floats(size = 3)

    e3 = 0. 0. 1.
      .type = floats(size = 3)

    relative_to_static_orientation = True
      .type = bool
  }
''')

help_message = '''

Generate plots of the scan-varying crystal orientation and unit cell from the
input refined_experiments.json which must contain a scan-varying crystal model.

Examples::

  dials.plot_scan_varying_crystal refined_experiments.json

'''

class Script(object):
  '''Class to run script.'''

  def __init__(self):
    '''Setup the script.'''
    from dials.util.options import OptionParser
    import libtbx.load_env

    usage = "usage: %s [options] experiments.json" % libtbx.env.dispatcher_name
    self.parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

  def run(self):
    '''Run the script.'''
    from dials.util.options import flatten_experiments
    from scitbx import matrix

    params, options = self.parser.parse_args()
    if len(params.input.experiments) == 0:
      self.parser.print_help()
      raise Sorry("No experiments found in the input")
    experiments = flatten_experiments(params.input.experiments)

    # Determine output path
    self._directory = os.path.join(params.output.directory,
      "scan-varying_crystal")
    self._directory = os.path.abspath(self._directory)
    ensure_directory(self._directory)
    self._format = "." + params.output.format

    self._debug = params.output.debug

    # Decomposition axes
    self._e1 = params.orientation_decomposition.e1
    self._e2 = params.orientation_decomposition.e2
    self._e3 = params.orientation_decomposition.e3

    # cell plot
    dat = []
    for iexp, exp in enumerate(experiments):

      crystal = exp.crystal
      scan = exp.scan

      if crystal.num_scan_points == 0:
        print("Ignoring scan-static crystal")
        continue

      scan_pts = range(crystal.num_scan_points)
      cells = [crystal.get_unit_cell_at_scan_point(t) for t in scan_pts]
      cell_params = [e.parameters() for e in cells]
      a,b,c,aa,bb,cc = zip(*cell_params)
      phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
      vol = [e.volume() for e in cells]
      cell_dat = {'phi':phi,
                  'a':a,
                  'b':b,
                  'c':c,
                  'alpha':aa,
                  'beta':bb,
                  'gamma':cc,
                  'volume':vol}
      if self._debug:
        print("Crystal in Experiment {0}".format(iexp))
        print("Phi\ta\tb\tc\talpha\tbeta\tgamma\tVolume")
        msg = "{0}\t{1}\t{2}\t{3}\t{4}\t{5}\t{6}\t{7}"
        line_dat = zip(phi, a, b, c, aa, bb, cc, vol)
        for line in line_dat:
          print(msg.format(*line))
      dat.append(cell_dat)
    self.plot_cell(dat)

    # orientation plot
    dat = []
    for iexp, exp in enumerate(experiments):

      crystal = exp.crystal
      scan = exp.scan

      if crystal.num_scan_points == 0:
        print("Ignoring scan-static crystal")
        continue

      scan_pts = range(crystal.num_scan_points)
      phi = [scan.get_angle_from_array_index(t) for t in scan_pts]
      Umats = [matrix.sqr(crystal.get_U_at_scan_point(t)) for t in scan_pts]
      if params.orientation_decomposition.relative_to_static_orientation:
        # factor out static U
        Uinv = matrix.sqr(crystal.get_U()).inverse()
        Umats = [U*Uinv for U in Umats]
      # NB e3 and e1 definitions for the crystal are swapped compared
      # with those used inside the solve_r3_rotation_for_angles_given_axes
      # method
      angles = [solve_r3_rotation_for_angles_given_axes(U,
        self._e3, self._e2, self._e1, deg=True) for U in Umats]
      phi3, phi2, phi1 = zip(*angles)
      angle_dat = {'phi':phi,
                   'phi3':phi3,
                   'phi2':phi2,
                   'phi1':phi1}
      if self._debug:
        print("Crystal in Experiment {0}".format(iexp))
        print("Image\tphi3\tphi2\tphi1")
        msg = "{0}\t{1}\t{2}\t{3}"
        line_dat = zip(phi, phi3, phi2, phi1)
        for line in line_dat:
          print(msg.format(*line))
      dat.append(angle_dat)
    self.plot_orientation(dat)

  def plot_cell(self, dat):
    try:
      import matplotlib.pyplot as plt
      import matplotlib.gridspec as gridspec
    except ImportError as e:
      print("matplotlib modules not available", e)
      return None

    from math import floor, ceil
    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(4, 2, wspace=0.4, hspace=0.6)

    ax = plt.subplot(gs[0, 0])
    ax.ticklabel_format(useOffset=False)
    for cell in dat:
      plt.plot(cell['phi'], cell['a'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('a')

    ax = plt.subplot(gs[0, 1])
    ax.ticklabel_format(useOffset=False)
    ymin, ymax = 0.0, 0.0
    for cell in dat:
      plt.plot(cell['phi'], cell['alpha'])
      # choose the widest y range
      ymin = max(ymin, floor(min(cell['alpha']) - 0.05))
      ymax = max(ymax, ceil(max(cell['alpha']) + 0.05))
      plt.axis(ymin=ymin, ymax=ymax)
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\alpha$')

    ax = plt.subplot(gs[1, 0])
    ax.ticklabel_format(useOffset=False)
    for cell in dat:
      plt.plot(cell['phi'], cell['b'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('b')

    ax = plt.subplot(gs[1, 1])
    ax.ticklabel_format(useOffset=False)
    ymin, ymax = 0.0, 0.0
    for cell in dat:
      plt.plot(cell['phi'], cell['beta'])
      # choose the widest y range
      ymin = max(ymin, floor(min(cell['beta']) - 0.05))
      ymax = max(ymax, ceil(max(cell['beta']) + 0.05))
      plt.axis(ymin=ymin, ymax=ymax)
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\beta$')

    ax = plt.subplot(gs[2, 0])
    ax.ticklabel_format(useOffset=False)
    for cell in dat:
      plt.plot(cell['phi'], cell['c'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('c')

    ax = plt.subplot(gs[2, 1])
    ax.ticklabel_format(useOffset=False)
    ymin, ymax = 0.0, 0.0
    for cell in dat:
      plt.plot(cell['phi'], cell['gamma'])
      # choose the widest y range
      ymin = max(ymin, floor(min(cell['gamma']) - 0.05))
      ymax = max(ymax, ceil(max(cell['gamma']) + 0.05))
      plt.axis(ymin=ymin, ymax=ymax)
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\gamma$')

    ax = plt.subplot2grid((4,2), (3, 0), colspan=2)
    ax.ticklabel_format(useOffset=False)
    for cell in dat:
      plt.plot(cell['phi'], cell['volume'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'volume $\left(\AA^3\right)$')
    plt.title('Cell volume')

    basename = os.path.join(self._directory, "unit_cell")
    fullname = basename + self._format
    print("Saving unit cell plot to {0}".format(fullname))
    plt.savefig(fullname)

  def plot_orientation(self, dat):
    try:
      import matplotlib.pyplot as plt
      import matplotlib.gridspec as gridspec
    except ImportError as e:
      print("matplotlib modules not available", e)
      return None

    from math import floor, ceil
    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(3, 1, wspace=0.4, hspace=0.6)

    ax = plt.subplot(gs[0, 0])
    ax.ticklabel_format(useOffset=False)
    for ori in dat:
      plt.plot(ori['phi'], ori['phi1'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_1$')

    ax = plt.subplot(gs[1, 0])
    ax.ticklabel_format(useOffset=False)
    for ori in dat:
      plt.plot(ori['phi'], ori['phi2'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_2$')

    ax = plt.subplot(gs[2, 0])
    ax.ticklabel_format(useOffset=False)
    for ori in dat:
      plt.plot(ori['phi'], ori['phi3'])
    plt.xlabel(r'rotation angle $\left(^\circ\right)$')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_3$')

    basename = os.path.join(self._directory, "orientation")
    fullname = basename + self._format
    print("Saving orientation plot to {0}".format(fullname))
    plt.savefig(fullname)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
