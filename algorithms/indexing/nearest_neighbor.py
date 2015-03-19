from __future__ import division
import math

class neighbor_analysis(object):
  def __init__(self, rs_vectors, tolerance=1.5, max_height_fraction=0.25,
               percentile=0.05):
    self.tolerance = tolerance # Margin of error for max unit cell estimate
    from scitbx.array_family import flex
    NEAR = 10
    self.NNBIN = 5 # target number of neighbors per histogram bin

    # nearest neighbor analysis
    from annlib_ext import AnnAdaptor
    query = flex.double()
    for spot in rs_vectors: # spots, in reciprocal space xyz
      query.append(spot[0])
      query.append(spot[1])
      query.append(spot[2])

    assert len(rs_vectors)>NEAR, (
      "Too few spots (%d) for nearest neighbour analysis." %len(rs_vectors))

    IS_adapt = AnnAdaptor(data=query,dim=3,k=1)
    IS_adapt.query(query)

    direct = flex.double()
    for i in xrange(len(rs_vectors)):
      direct.append(1.0/math.sqrt(IS_adapt.distances[i]))
    # reject top 1% of longest distances to hopefully get rid of any outliers
    direct = direct.select(
      flex.sort_permutation(direct)[:int(math.floor(0.99*len(direct)))])

    # determine the most probable nearest neighbor distance (direct space)
    hst = flex.histogram(direct, n_slots=int(len(direct)/self.NNBIN))
    centers = hst.slot_centers()
    islot = hst.slots()
    highest_bin_height = flex.max(islot)
    most_probable_neighbor = centers[list(islot).index(highest_bin_height)]

    if False:  # to print out the histogramming analysis
      smin, smax = flex.min(direct), flex.max(direct)
      stats = flex.mean_and_variance(direct)
      import sys
      out = sys.stdout
      print >> out, "     range:     %6.2f - %.2f" % (smin, smax)
      print >> out, "     mean:      %6.2f +/- %6.2f on N = %d" % (
        stats.mean(), stats.unweighted_sample_standard_deviation(), direct.size())
      hst.show(f=out, prefix="    ", format_cutoffs="%6.2f")
      print >> out, ""

    # choose a max cell based on bins above a given fraction of the highest bin height
    # given multiple
    isel = (islot.as_double() > (
      max_height_fraction * highest_bin_height)).iselection()
    self.max_cell = (
      self.tolerance * centers[int(flex.max(isel.as_double()))]+0.5*hst.slot_width())

    # determine the 5th-percentile direct-space distance
    perm = flex.sort_permutation(direct, reverse=True)
    self.percentile = direct[perm[int(percentile * len(rs_vectors))]]

    #MAXTOL = 1.5 # Margin of error for max unit cell estimate
    #self.max_cell = max( MAXTOL * most_probable_neighbor,
                         #MAXTOL * self.percentile)

    if False:
      self.plot(direct)

  def plot(self,val):
    import numpy as np

    hist,bins = np.histogram(val,bins=int(len(val)/self.NNBIN))
    width = 0.7*(bins[1]-bins[0])
    center = (bins[:-1]+bins[1:])/2
    import matplotlib.pyplot as plt
    plt.bar(center, hist, align="center", width=width)
    ymin, ymax = plt.ylim()
    plt.vlines(self.percentile, ymin, ymax, colors='r')
    plt.vlines(self.max_cell/self.tolerance, ymin, ymax, colors='g')
    plt.show()
