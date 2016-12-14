from __future__ import division
import math

class neighbor_analysis(object):
  def __init__(self, reflections, step_size=45, tolerance=1.5,
               max_height_fraction=0.25, percentile=0.05):
    self.tolerance = tolerance # Margin of error for max unit cell estimate
    from scitbx.array_family import flex
    NEAR = 10
    self.NNBIN = 5 # target number of neighbors per histogram bin

    direct = flex.double()

    if 'entering' in reflections:
      entering_flags = reflections['entering']
    else:
      entering_flags = flex.bool(reflections.size(), True)
    rs_vectors = reflections['rlp']
    phi_deg = reflections['xyzobs.mm.value'].parts()[2] * (180/math.pi)

    d_spacings = flex.double()
    # nearest neighbor analysis
    from annlib_ext import AnnAdaptor
    for imageset_id in range(flex.max(reflections['imageset_id'])+1):
      sel = reflections['imageset_id'] == imageset_id
      phi_min = flex.min(phi_deg.select(sel))
      phi_max = flex.max(phi_deg.select(sel))
      d_phi = phi_max - phi_min
      n_steps = max(int(math.ceil(d_phi / step_size)), 1)

      for n in range(n_steps):
        sel &= (phi_deg >= (phi_min+n*step_size)) & (phi_deg < (phi_min+(n+1)*step_size))

        for entering in (True, False):
          sel  &= entering_flags == entering
          if sel.count(True) == 0:
            continue

          query = flex.double()
          query.extend(rs_vectors.select(sel).as_double())

          if query.size() == 0:
            continue

          IS_adapt = AnnAdaptor(data=query,dim=3,k=1)
          IS_adapt.query(query)

          direct.extend(1/flex.sqrt(IS_adapt.distances))
          d_spacings.extend(1/rs_vectors.norms())

    assert len(direct)>NEAR, (
      "Too few spots (%d) for nearest neighbour analysis." %len(direct))

    perm = flex.sort_permutation(direct)
    direct = direct.select(perm)
    d_spacings = d_spacings.select(perm)

    # reject top 1% of longest distances to hopefully get rid of any outliers
    n = int(math.floor(0.99*len(direct)))
    direct = direct[:n]
    d_spacings = d_spacings[:n]

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
    self.percentile = direct[perm[int(percentile * len(direct))]]

    #MAXTOL = 1.5 # Margin of error for max unit cell estimate
    #self.max_cell = max(MAXTOL * most_probable_neighbor,
                         #MAXTOL * self.percentile)

    self.reciprocal_lattice_vectors = rs_vectors
    self.d_spacings = d_spacings
    self.direct = direct
    self.histogram = hst

  def plot_histogram(self, filename='nn_hist.png', figsize=(12,8)):
    import matplotlib.pyplot as plt
    hist = self.histogram
    fig = plt.figure(figsize=figsize)
    plt.bar(hist.slot_centers(), hist.slots(), align="center",
            width=hist.slot_width(), color='black', edgecolor=None)
    ymin, ymax = plt.ylim()
    plt.vlines(self.max_cell/self.tolerance, ymin, ymax, colors='g', label='estimated max cell')
    plt.xlabel('Direct space distance (A)')
    plt.ylabel('Frequency')
    plt.legend(loc='best')
    plt.savefig(filename)
    plt.clf()
