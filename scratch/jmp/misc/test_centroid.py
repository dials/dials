



from dials.model.serialize import load
spots2 = load.reflections('strong2.pickle')
print spots2[52776]
#for i, s in enumerate(spots2):
#    x, y, z = s.centroid_position
#    if z >= 357 and z < 372:
#        if x >= 1240 and x < 1250:
#            print i, x, y, z
#    if i == 49807:
#        print s

spots = load.reflections('strong.pickle')

#dx, dy, dz = [], [], []
#for s1, s2 in zip(spots, spots2):
#    dx.append(s1.centroid_variance[0] - s2.centroid_variance[0])
#    dx.append(s1.centroid_variance[1] - s2.centroid_variance[1])
#    dx.append(s1.centroid_variance[2] - s2.centroid_variance[2])

#from matplotlib import pylab
#pylab.plot(dx)
#pylab.show()

#for s1, s2 in zip(spots, spots2):
#    print s1.centroid_variance, s2.centroid_variance

#from dials.algorithms.centroid import ComputeCentroid

#compute_centroid = ComputeCentroid()
#compute_centroid(spots)

def hist_outline(hist):
  from scitbx.array_family import flex
  step_size = hist.slot_width()
  half_step_size = 0.5 * step_size
  n_slots = len(hist.slots())

  bins = flex.double(n_slots * 2 + 2, 0)
  data = flex.double(n_slots * 2 + 2, 0)
  for i in range(n_slots):
    bins[2 * i + 1] = hist.slot_centers()[i] - half_step_size
    bins[2 * i + 2] = hist.slot_centers()[i] + half_step_size
    data[2 * i + 1] = hist.slots()[i]
    data[2 * i + 2] = hist.slots()[i]

  bins[0] = bins[1] - step_size
  bins[-1] = bins[-2] + step_size
  data[0] = 0
  data[-1] = 0

  return (bins, data)

def plot_centroid_weights_histograms(reflections, n_slots=50):
  from matplotlib import pyplot
  from scitbx.array_family import flex
  variances = flex.vec3_double([r.centroid_variance for r in reflections])
  vx, vy, vz = variances.parts()
  idx = (vx > 0).__and__(vy > 0).__and__(vz > 0)
  vx = vx.select(idx)
  vy = vy.select(idx)
  vz = vz.select(idx)
  wx = 1/vx
  wy = 1/vy
  wz = 1/vz
  wx = flex.log(wx)
  wy = flex.log(wy)
  wz = flex.log(wz)
  hx = flex.histogram(wx, n_slots=n_slots)
  hy = flex.histogram(wy, n_slots=n_slots)
  hz = flex.histogram(wz, n_slots=n_slots)
  fig = pyplot.figure()

  idx2 = flex.max_index(wx)
  idx3 = flex.int(range(len(reflections))).select(idx)[idx2]
  print reflections[idx3]
  return

  #outliers = reflections.select(wx > 50)
  #for refl in outliers:
    #print refl

  for i, h in enumerate([hx, hy, hz]):
    ax = fig.add_subplot(311+i)

    slots = h.slots().as_double()
    bins, data = hist_outline(h)
    log_scale = True
    if log_scale:
      data.set_selected(data == 0, 0.1) # otherwise lines don't get drawn when we have some empty bins
      ax.set_yscale("log")
    ax.plot(bins, data, '-k', linewidth=2)
    #pyplot.suptitle(title)
    data_min = min([slot.low_cutoff for slot in h.slot_infos() if slot.n > 0])
    data_max = max([slot.low_cutoff for slot in h.slot_infos() if slot.n > 0])
    ax.set_xlim(data_min, data_max+h.slot_width())
  pyplot.show()


plot_centroid_weights_histograms(spots)
