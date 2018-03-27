#!/usr/bin/env python
#
# plot_centroid_analysis.py
#
#  Copyright (C) (2017) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David G. Waterman
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.

# LIBTBX_SET_DISPATCHER_NAME dev.dials.plot_centroid_analysis

from __future__ import absolute_import, division
from __future__ import print_function
from libtbx.utils import Sorry
from libtbx.phil import parse
from dials.array_family import flex
from dials.algorithms.refinement.analysis.centroid_analysis import CentroidAnalyser

phil_scope = parse('''

  plot {
    mark_periods = 54.0 36.0 18.0
      .type = floats(value_min=1)
      .help = "Draw vertical lines at frequencies corresponding to these periods"
              "in degrees"
  }

  output {
    format = *png pdf
      .type = choice
  }
''')

help_message = '''

Generate plots and tables of centroid analysis of residuals. This centroid
analysis is used by dials.refine in order to set a suitable block width for
outlier rejection and to set an appropriate interval_width_degrees parameter
for scan-varying refinement.

Examples::

  dev.dials.plot_centroid_analysis reflections.pickle

'''

def save_plots(params, raw, smoothed, suffix=''):
  """Create plots for the centroid analysis results for a single experiment.
  Overlay raw and smoothed periodograms"""

  # draw vertical lines at frequencies corresponding to periods of 54, 36 and 18
  # degrees
  if params.plot.mark_periods is not None:
    vlines = [1./e for e in params.plot.mark_periods]
  else:
    vlines = []

  import matplotlib
  matplotlib.use('Agg')
  import matplotlib.pyplot as plt

  nblocks = raw['nblocks']
  block_size = raw['block_size']
  phistart = raw['phi_range'][0]
  block_centres = block_size * flex.double_range(nblocks) + phistart + block_size/2.0

  # X residuals plot
  plt.figure(1)
  plt.subplot(211)
  plt.plot(block_centres, 1000. * raw['av_x_resid_per_block'])
  plt.xlabel('phi (degrees)')
  plt.ylabel('x residuals per block (microns)')

  # X periodogram
  plt.subplot(212)
  for dat in [raw, smoothed]: # overlay raw and smoothed periodogram plots
    px = dat['x_periodogram']
    if px is None: continue
    sample_freq = 1./dat['block_size']
    freq = px.freq * sample_freq
    line, = plt.semilogy(freq, px.spec)
  for vline in vlines:
    plt.axvline(x=vline, color='r')
  x_interval = smoothed['x_interval']
  if x_interval:
    x_freq = 1./x_interval
    plt.axvline(x=x_freq, color='c')
    plt.text(0.2, 0.9, 'interval width: {0:.3f}'.format(x_interval),
             transform=line.axes.transAxes, color='c')
  plt.xlabel('frequency')
  plt.ylabel('spectrum')

  # write them out
  fname = 'x-residual-analysis' + suffix + '.' + params.output.format
  print("Saving {0}".format(fname))
  plt.savefig(fname)

  # Y residuals plot
  plt.figure(2)
  plt.subplot(211)
  plt.plot(block_centres, 1000. * raw['av_y_resid_per_block'])
  plt.xlabel('phi (degrees)')
  plt.ylabel('y residuals per block (microns)')

  # Y periodogram
  plt.subplot(212)
  for dat in [raw, smoothed]: # overlay raw and smoothed periodogram plots
    py = dat['y_periodogram']
    if py is None: continue
    sample_freq = 1./dat['block_size']
    freq = py.freq * sample_freq
    line, = plt.semilogy(freq, py.spec)
  for vline in vlines:
    plt.axvline(x=vline, color='r')
  y_interval = smoothed['y_interval']
  if y_interval:
    y_freq = 1./y_interval
    plt.axvline(x=y_freq, color='c')
    plt.text(0.2, 0.9, 'interval width: {0:.3f}'.format(y_interval),
             transform=line.axes.transAxes, color='c')
  plt.xlabel('frequency')
  plt.ylabel('spectrum')

  # write them out
  fname = 'y-residual-analysis' + suffix + '.' + params.output.format
  print("Saving {0}".format(fname))
  plt.savefig(fname)

  # phi residuals plot
  plt.figure(3)
  plt.subplot(211)
  plt.plot(block_centres, 1000. * raw['av_phi_resid_per_block'])
  plt.xlabel('phi (degrees)')
  plt.ylabel('phi residuals per block (mrad)')

  # phi periodogram
  plt.subplot(212)
  for dat in [raw, smoothed]: # overlay raw and smoothed periodogram plots
    pz = dat['phi_periodogram']
    if pz is None: continue
    sample_freq = 1./dat['block_size']
    freq = pz.freq * sample_freq
    line, = plt.semilogy(freq, pz.spec)
  for vline in vlines:
    plt.axvline(x=vline, color='r')
  phi_interval = smoothed['phi_interval']
  if phi_interval:
    phi_freq = 1./phi_interval
    plt.axvline(x=phi_freq, color='c')
    plt.text(0.2, 0.9, 'interval width: {0:.3f}'.format(phi_interval),
             transform=line.axes.transAxes, color='c')
  plt.xlabel('frequency')
  plt.ylabel('spectrum')

  # write them out
  fname = 'phi-residual-analysis' + suffix + '.' + params.output.format
  print("Saving {0}".format(fname))
  plt.savefig(fname)

  return

def run(args):

  import libtbx.load_env
  usage = """\
%s reflections.pickle [options]""" %libtbx.env.dispatcher_name
  from dials.util.options import OptionParser
  from dials.util.options import flatten_reflections
  from libtbx.utils import Sorry
  parser = OptionParser(
    usage=usage,
    phil=phil_scope,
    read_reflections=True,
    epilog=help_message)

  params, options = parser.parse_args(show_diff_phil=True)
  reflections = flatten_reflections(params.input.reflections)

  if len(reflections) != 1:
    parser.print_help()
    raise Sorry("Please provide a single file of reflections")

  refs = reflections[0]

  # raw periodograms
  ca = CentroidAnalyser(refs)
  results_r = ca(spans=None)

  # smoothed periodograms
  ca = CentroidAnalyser(refs)
  results_s = ca()

  if len(results_r) == 1:
    save_plots(params, results_r[0], results_s[0])
  else:
    for i, (r_r, r_s) in enumerate(zip(results_r, results_s)):
      suffix = '_exp_{0}'.format(i)
      save_plots(params, r_r, r_s, suffix=suffix)

  # TODO: print tables of data from the analysis

if __name__ == '__main__':
  import sys
  run(sys.argv[1:])
