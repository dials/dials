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


class ScriptRunner(object):
  '''Class to run script.'''

  def __init__(self, crystals):
    '''Setup the script.'''

    # Filename data
    self.crystals = crystals

  def __call__(self):
    '''Run the script.'''

    for icrystal, crystal in enumerate(self.crystals):

      if crystal.num_scan_points == 0:
        print "Ignoring scan-static crystal"
        continue

      scan_pts = range(crystal.num_scan_points)
      cells = [crystal.get_unit_cell_at_scan_point(t) \
               for t in scan_pts]
      dat = [(t,) + e.parameters() + (e.volume(),) \
             for (t, e) in zip(scan_pts, cells)]
      self.plot(dat)

      print "Image\ta\tb\tc\talpha\tbeta\tgamma\tVolume"
      msg = "\t".join(["%.3f"] * 8)
      for line in dat:
        print msg % line

    print "TODO: misset angles around user-supplied axes"

  def plot(self, dat):
    try:
      import matplotlib.pyplot as plt
      import matplotlib.gridspec as gridspec
    except ImportError as e:
      print "matplotlib modules not available", e
      return None

    from math import floor, ceil
    image, a, b, c, alpha, beta, gamma, volume = zip(*dat)

    gs = gridspec.GridSpec(4, 2, wspace=0.4, hspace=0.6)
    plt.subplot(gs[0, 0])
    plt.plot(image, a)
    plt.xlabel('Image')
    plt.ylabel('$\AA$')
    plt.title('a')

    plt.subplot(gs[0, 1])
    plt.plot(image, alpha)
    plt.axis(ymin=floor(min(alpha)), ymax=ceil(max(alpha)))
    plt.xlabel('Image')
    plt.ylabel(r'$\AA$')
    plt.title(r'$\alpha$')

    plt.subplot(gs[1, 0])
    plt.plot(image, b)
    plt.xlabel('Image')
    plt.ylabel(r'$\AA$')
    plt.title('b')

    plt.subplot(gs[1, 1])
    plt.plot(image, beta)
    plt.axis(ymin=floor(min(beta)), ymax=ceil(max(beta)))
    plt.xlabel('Image')
    plt.ylabel(r'$\AA$')
    plt.title(r'$\beta$')

    plt.subplot(gs[2, 0])
    plt.plot(image, c)
    plt.xlabel('Image')
    plt.ylabel(r'$\AA$')
    plt.title('c')

    plt.subplot(gs[2, 1])
    plt.plot(image, gamma)
    plt.axis(ymin=floor(min(gamma)), ymax=ceil(max(gamma)))
    plt.xlabel('Image')
    plt.ylabel(r'$\AA$')
    plt.title(r'$\gamma$')

    plt.subplot2grid((4,2), (3, 0), colspan=2)
    plt.plot(image, volume)
    plt.xlabel('Image')
    plt.ylabel(r'$\AA^3$')
    plt.title('Cell volume')

    plt.show()

if __name__ == '__main__':
  import sys
  from dials.util.command_line import Importer
  args = sys.argv[1:]
  importer = Importer(args, check_format=False)
  try:
    crystals = importer.experiments.crystals()
  except AttributeError:
    print "No crystals found in the input"
    raise

  assert len(importer.unhandled_arguments) == 0

  runner = ScriptRunner(
      crystals=crystals)

  # Run the script
  runner()
