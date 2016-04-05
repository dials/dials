#!/usr/bin/env cctbx.python

#
#  Copyright (C) (2016) STFC Rutherford Appleton Laboratory, UK.
#
#  Author: David Waterman.
#
#  This code is distributed under the BSD license, a copy of which is
#  included in the root directory of this package.
#

# LIBTBX_SET_DISPATCHER_NAME dev.dials.bijvoet_beam_search
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export PHENIX_GUI_ENVIRONMENT=1
# LIBTBX_PRE_DISPATCHER_INCLUDE_SH export BOOST_ADAPTBX_FPE_DEFAULT=1


"""
Proof of concept for an algorithm to search for the beam centre (assuming a
single flat panel detector) using the intersection of lines bisecting vectors
joining supposed symmetry-related pairs of reflections in detector space.

"""

from __future__ import division

from dxtbx.model.experiment.experiment_list import Experiment
from dials.util.options import flatten_reflections
from scitbx.array_family import flex
from scitbx import matrix
from libtbx import phil, Auto
from libtbx.utils import Sorry
import matplotlib
# Offline backend
#matplotlib.use("Agg")
matplotlib.rc('font', family='serif')
matplotlib.rc('font', serif='Times New Roman')
from matplotlib import pyplot as plt

help_message = '''

Plots lines between reflections with similar intensities (which are possibly
Bijvoet mates) and the bisectors of these lines, which converge close to
the beam centre if the reflections are true symmetry mates.

Example::

  dev.dials.bijvoet_beam_search datablock.json strong.pickle

'''

phil_scope = phil.parse('''
num_reflections = 300
  .type = int(value_min = 2)
  .help = "The number of strongest reflections to calculate potential Bijvoet"
          "pairings for."
pair_cutoff
{
  frac_intensity = 0.001
    .type = float
    .help = "Fractional intensity difference below which a potential Bijvoet"
            "pair may be constructed. f = 0.5 * |I1 - I2| / (I1 + I2). If"
            "f < pair_cutoff, assume a Bijvoet pairing. In future, better to"
            "use probabilities"
  distance = 10
    .type = float
    .help = "Distance cutoff in millimetres below which two centroids will"
            "not be considered a Bijvoet pair"
}
beam_centre_convention = *mosflm dials
  .type = choice
  .help = "Convention for the detector mm coordinate system. By default"
          "set to mosflm to easily set the mosflm_beam_centre parameter"
          "of dials.import"
plot
{
  bijvoet_connectors = False
    .type = bool
    .help = "Whether to plot the positions of the supposed Bijvoet mates and"
            "the connectors between them"
}
''')

class PutativeBijvoetPair(object):

  def __init__(self, r1, r2):

    x1, y1 = r1['detX'], r1['detY']
    x2, y2 = r2['detX'], r2['detY']

    self.x1 = x1
    self.y1 = y1
    self.x2 = x2
    self.y2 = y2

    self.I1 = r1['intensity.sum.value']
    self.I2 = r2['intensity.sum.value']

    pt1 = matrix.col((x1, y1))
    pt2 = matrix.col((x2, y2))

    connector = (pt2 - pt1)
    self.distance = connector.length()

    midpoint = pt1 + 0.5 * connector

    # unit vector along the bisector
    R = matrix.sqr((0, -1, 1, 0))
    ubi = (R * connector).normalize()
    # A and B are points along the bisector line hard coded at 1 metre away
    # from the midpoint
    A = midpoint + 1000 * ubi
    B = midpoint - 1000 * ubi

    self.bisector_x1 = A[0]
    self.bisector_y1 = A[1]
    self.bisector_x2 = B[0]
    self.bisector_y2 = B[1]

    return

class PutativeBijvoetPairFactory(object):

  def __init__(self, intensity_cutoff, distance_cutoff):

    self.intensity_cutoff = intensity_cutoff
    self.distance_cutoff = distance_cutoff
    self.last_refusal_reason = None

  def __call__(self, r1, r2):

    pair = PutativeBijvoetPair(r1, r2)
    if pair.distance < self.distance_cutoff:
      self.last_refusal_reason = 'distance'
      return None

    fracIdiff = 0.5 * abs(pair.I1 - pair.I2) / (pair.I1 + pair.I2)
    if fracIdiff >= self.intensity_cutoff:
      self.last_refusal_reason = 'intensity'
      return None

    return pair


class Script(object):

  def __init__(self):
    '''Check script input'''

    import libtbx.load_env
    from dials.util.options import OptionParser

    # The script usage
    usage = ("usage: {0} [options] [param.phil] experiments1.json "
             "experiments2.json").format(libtbx.env.dispatcher_name)

    parser = OptionParser(
      usage=usage,
      phil=phil_scope,
      read_reflections=True,
      read_datablocks=True,
      read_experiments=True,
      check_format=False,
      epilog=help_message)

    params, options = parser.parse_args(show_diff_phil=True)

    warnmsg = ("WARNING: The {0} contains more than one "
               "{1}. Only the first will be considered.")
    if params.input.datablock:
      if params.input.experiments:
        raise Sorry("Please provide either a datablock or an experiment list, "
                    "not a combination of the two")
      if len(params.input.datablock) != 1:
        raise Sorry("Please provide a single datablock or experiment list as input")

      detector = params.input.datablock[0].data[0].unique_detectors()
      if len(detector) > 1:
        print warnmsg.format("datablock", "detector")
      detector = detector[0]
      beam = params.input.datablock[0].data[0].unique_beams()
      if len(beam) > 1:
        print warnmsg.format("datablock", "beam")
      beam = beam[0]
      scan = params.input.datablock[0].data[0].unique_scans()
      if len(scan) > 1:
        print warnmsg.format("datablock", "scan")
      scan = scan[0]

    else:
      if len(params.input.experiments) != 1:
        raise Sorry("Please provide a single datablock or experiment list as input")

      detector = params.input.experiments[0].data.detectors()
      if len(detector) > 1:
        print warnmsg.format("experiment list", "detector")
      detector = detector[0]

      beam = params.input.experiments[0].data.beams()
      if len(beam) > 1:
        print warnmsg.format("experiment list", "beam")
      beam = beam[0]
      scan = params.input.experiments[0].data.scans()
      if len(scan) > 1:
        print warnmsg.format("experiment list", "scan")
      scan = scan[0]

    self.detector = detector
    self.beam = beam
    self.scan = scan
    self.params = params

    reflections = flatten_reflections(params.input.reflections)
    if len(reflections) > 1:
      raise Sorry("Please provide a single reflection table as input")
    self.reflections = reflections[0]
    if not self.reflections.has_key('intensity.sum.value'):
      raise Sorry('The provided reflection table does not have a column '
                  'named "intensity.sum.value". Cannot continue.')
    self.reflections.sort(name='intensity.sum.value', reverse=True)

    from dials.algorithms.indexing import indexer
    self.reflections = indexer.indexer_base.map_spots_pixel_to_mm_rad(
      self.reflections, self.detector, self.scan)

    if len(self.detector) != 1:
      raise Sorry('The detector has multiple panels, which is not currently '
                  'supported by this program.')
    self.panel = detector[0]

    # For convenience, set centroid X and Y in their own columns of the
    # reflection table, taking into account the coordinate system convention.
    # If using the mosflm coordinate system this has X and Y reversed wrt DIALS.
    X, Y, phi = self.reflections['xyzobs.mm.value'].parts()
    if params.beam_centre_convention == 'mosflm':
      (X, Y) = (Y, X)
    else: assert self.params.beam_centre_convention == 'dials'
    self.reflections['detX'] = X
    self.reflections['detY'] = Y

    self.pair_factory = PutativeBijvoetPairFactory(
      intensity_cutoff=params.pair_cutoff.frac_intensity,
      distance_cutoff=params.pair_cutoff.distance)

    return

  def __call__(self):

    # Take the top most intense reflections
    top = self.reflections[0:self.params.num_reflections]

    # Set up plot
    fig=plt.figure()
    plt.xlabel("X (mm)")
    plt.ylabel("Y (mm)")
    ax=fig.add_subplot(111)
    ax.set_aspect(1)
    ax.set_title(('Bisectors of lines joining potential Bijvoet mates.\n'
                  'X and Y using the {0} convention').format(
                   self.params.beam_centre_convention.upper()))

    # Start with the most intense spot, consider whether it is a symmetry
    # mate of the next spot, and so on down the list until the it no longer
    # seems 'likely' that these are true pairs (XXX This part needs more work
    # really). Then move on to the second most intense spot and work downwards
    # as before. Each time, keep a PutativeBijvoetPair, containing data for the
    # reflection positions, the line between them and the bisector line.
    pairs = []
    for i, r1 in enumerate(top):
      for r2 in top[(i+1):]:
        pair = self.pair_factory(r1, r2)
        if pair is None:
          # as soon as intensity is a reason for refusing to construct a
          # Bijvoet pair, give up as later reflections are even more different
          if self.pair_factory.last_refusal_reason == 'intensity':
            break
        else: pairs.append(pair)

    for pair in pairs:
      if self.params.plot.bijvoet_connectors:
        ax.plot([pair.x1, pair.x2],[pair.y1, pair.y2],"ro")
        ax.plot([pair.x1, pair.x2],[pair.y1, pair.y2],"r--")

      ax.plot([pair.bisector_x1, pair.bisector_x2],
              [pair.bisector_y1, pair.bisector_y2],
              'g-', alpha=0.1, linewidth=5)

    # get current beam centre and detector size
    beam_centre = self.panel.get_ray_intersection(self.beam.get_s0())
    xlim, ylim = self.panel.get_image_size_mm()
    if self.params.beam_centre_convention == 'mosflm':
      beam_centre = beam_centre[1], beam_centre[0]
      xlim, ylim = ylim, xlim

    # plot the current beam centre
    ax.plot(beam_centre[0], beam_centre[1], 'bo')

    # set a sensible plot window field of view
    #ax.set_xlim(beam_centre[0] - xlim / 20, beam_centre[0] + xlim / 20)
    #ax.set_ylim(beam_centre[1] - ylim / 20, beam_centre[1] + ylim / 20)
    ax.set_xlim(beam_centre[0] - 10, beam_centre[0] + 10)
    ax.set_ylim(beam_centre[1] - 10, beam_centre[1] + 10)

    plt.show()
    plt.close()

if __name__ == "__main__":

  run = Script()
  run()
