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

from __future__ import division

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

class Script(object):
  '''Class to run script.'''

  def __init__(self):
    '''Setup the script.'''
    from dials.util.options import OptionParser

    self.parser = OptionParser(
      phil=phil_scope,
      read_experiments=True,
      check_format=False)

  def run(self):
    '''Run the script.'''
    from dials.util.options import flatten_experiments

    params, options = self.parser.parse_args()
    if len(params.input.experiments) == 0:
      self.parser.print_help()
      raise Sorry("No experiments found in the input")
    experiments = flatten_experiments(params.input.experiments)
    crystals = experiments.crystals()

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

    for icrystal, crystal in enumerate(crystals):

      icrystal_suffix = icrystal if len(crystals) > 1 else None

      if crystal.num_scan_points == 0:
        print "Ignoring scan-static crystal"
        continue

      # cell plot
      scan_pts = range(crystal.num_scan_points)
      cells = [crystal.get_unit_cell_at_scan_point(t) for t in scan_pts]
      dat = [(t,) + e.parameters() + (e.volume(),) \
             for (t, e) in zip(scan_pts, cells)]
      self.plot_cell(dat, icrystal_suffix)

      if self._debug:
        print "Crystal {0}".format(icrystal)
        print "Image\ta\tb\tc\talpha\tbeta\tgamma\tVolume"
        msg = "\t".join(["%.3f"] * 8)
        for line in dat:
          print msg % line

      # orientation plot
      Umats = [crystal.get_U_at_scan_point(t) for t in scan_pts]
      if params.orientation_decomposition.relative_to_static_orientation:
        # factor out static U
        Uinv = crystal.get_U().inverse()
        Umats = [U*Uinv for U in Umats]
      # NB e3 and e1 definitions for the crystal are swapped compared
      # with those used inside the solve_r3_rotation_for_angles_given_axes
      # method
      angles = [solve_r3_rotation_for_angles_given_axes(U,
        self._e3, self._e2, self._e1, deg=True) for U in Umats]
      dat = [(t,) + a for (t, a) in zip(scan_pts, angles)]
      self.plot_orientation(dat, icrystal_suffix)

      if self._debug:
        print "Crystal {0}".format(icrystal)
        print "Image\tphi3\tphi2\tphi1"
        msg = "\t".join(["%.6f"] * 4)
        for line in dat:
          print msg % line

  def plot_cell(self, dat, icrystal=None):
    try:
      import matplotlib.pyplot as plt
      import matplotlib.gridspec as gridspec
    except ImportError as e:
      print "matplotlib modules not available", e
      return None

    from math import floor, ceil
    image, a, b, c, alpha, beta, gamma, volume = zip(*dat)
    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(4, 2, wspace=0.4, hspace=0.6)

    ax = plt.subplot(gs[0, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, a)
    plt.xlabel('Image')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('a')

    ax = plt.subplot(gs[0, 1])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, alpha)
    plt.axis(ymin=floor(min(alpha) - 0.05), ymax=ceil(max(alpha) + 0.05))
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\alpha$')

    ax = plt.subplot(gs[1, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, b)
    plt.xlabel('Image')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('b')

    ax = plt.subplot(gs[1, 1])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, beta)
    plt.axis(ymin=floor(min(beta)- 0.05), ymax=ceil(max(beta) + 0.05))
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\beta$')

    ax = plt.subplot(gs[2, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, c)
    plt.xlabel('Image')
    plt.ylabel(r'length $\left(\AA\right)$')
    plt.title('c')

    ax = plt.subplot(gs[2, 1])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, gamma)
    plt.axis(ymin=floor(min(gamma) - 0.05), ymax=ceil(max(gamma) + 0.05))
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\gamma$')

    ax = plt.subplot2grid((4,2), (3, 0), colspan=2)
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, volume)
    plt.xlabel('Image')
    plt.ylabel(r'volume $\left(\AA^3\right)$')
    plt.title('Cell volume')

    basename = os.path.join(self._directory, "unit_cell")
    if icrystal is not None: basename += "_{0}".format(
      icrystal)
    plt.savefig(basename + self._format)

  def plot_orientation(self, dat, icrystal=None):
    try:
      import matplotlib.pyplot as plt
      import matplotlib.gridspec as gridspec
    except ImportError as e:
      print "matplotlib modules not available", e
      return None

    from math import floor, ceil
    image, phi3, phi2, phi1 = zip(*dat)
    fig = plt.figure(figsize=(13, 10))
    gs = gridspec.GridSpec(3, 1, wspace=0.4, hspace=0.6)

    ax = plt.subplot(gs[0, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, phi1)
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_1$')

    ax = plt.subplot(gs[1, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, phi2)
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_2$')

    ax = plt.subplot(gs[2, 0])
    ax.ticklabel_format(useOffset=False)
    plt.plot(image, phi3)
    plt.xlabel('Image')
    plt.ylabel(r'angle $\left(^\circ\right)$')
    plt.title(r'$\phi_3$')

    basename = os.path.join(self._directory, "orientation")
    if icrystal is not None: basename += "_{0}".format(
      icrystal)
    plt.savefig(basename + self._format)

if __name__ == '__main__':
  from dials.util import halraiser
  try:
    script = Script()
    script.run()
  except Exception as e:
    halraiser(e)
