
from __future__ import division

def fit(data, mask, model_data, model_mask, kernel):

  from dials.array_family import flex
  from dials.algorithms.image.filter import summed_area

  mask = mask & model_mask
  mask = mask.as_1d().as_int()
  mask.reshape(data.accessor())
  data = data * mask
  model_data = model_data * mask.as_double()

  summed_data = summed_area(data, kernel).as_1d()
  summed_mask = summed_area(mask, kernel).as_1d()
  summed_model = summed_area(model_data, kernel).as_1d()
  scaled_mask = (model_mask.as_1d() == True) & (summed_mask > 2) & (summed_model > 0)

  indices = flex.size_t(range(len(data))).select(scaled_mask)

  scaled_data = flex.double(len(data))
  scaled_data.set_selected(
    indices,
    summed_data.select(indices).as_double() / summed_model.select(indices))
  scaled_data.reshape(data.accessor())
  scaled_mask.reshape(data.accessor())

  #for j in range(scaled_data.all()[0]):
  #  last_scale = None
  #  for i in range(scaled_data.all()[1]):
  #    if scaled_mask[j,i] == False:
  #      if last_scale != None:
  #        scaled_data[j,i] = last_scale
  #        scaled_mask[j,i] = True
  #    else:
  #      last_scale = scaled_data[j,i]


  return scaled_data, scaled_mask





if __name__ == '__main__':

  from dials.util.options import OptionParser
  from dials.util.options import flatten_experiments
  from dials.util.options import flatten_reflections
  from dials.array_family import flex
  from collections import defaultdict
  from dials.algorithms.shoebox import MaskCode
  from libtbx.phil import parse
  import cPickle as pickle

  phil_scope = parse('''
    model = None
      .type = str
  ''')

  # Create the option parser
  parser = OptionParser(
    phil=phil_scope,
    read_experiments=True,
    read_reflections=True)

  # Parse the arguments
  params, options = parser.parse_args(show_diff_phil=True)

  assert(params.model is not None)
  model_data, model_mask = pickle.load(open(params.model))

  # Get the experiments and reflections
  experiments = flatten_experiments(params.input.experiments)
  reflections = flatten_reflections(params.input.reflections)

  # Get the imageset
  assert len(experiments) == 1
  assert len(reflections) == 1
  reflections = reflections[0]
  imageset = experiments[0].imageset

  print list(reflections.keys())
  reflections.unset_flags(flex.size_t(range(len(reflections))), reflections.flags.integrated_sum)
  reflections.unset_flags(flex.size_t(range(len(reflections))), reflections.flags.integrated_prf)
  del reflections['intensity.sum.value']
  del reflections['intensity.sum.variance']
  del reflections['background.mean']
  del reflections['intensity.prf.value']
  del reflections['intensity.prf.variance']
  del reflections['profile.correlation']
  del reflections['profile.rmsd']

  # Create the reflection lookup
  reflections['shoebox'] = flex.shoebox(
    reflections['panel'],
    reflections['bbox'])
  bbox = reflections['bbox']
  reflection_lookup = defaultdict(list)
  for i in range(len(bbox)):
    for j in range(bbox[i][4],bbox[i][5]):
      reflection_lookup[j].append(i)

  width, height = experiments[0].detector[0].get_image_size()

  # Loop through all images
  print "START"
  for frame in range(0, len(imageset)):

    # Get the subset of reflections on this image and compute the mask
    indices = flex.size_t(reflection_lookup[frame])
    x0, x1, y0, y1, z0, z1 = reflections['bbox'].parts()
    indices2 = flex.size_t(range(len(reflections))).select(z0 == frame)
    indices3 = flex.size_t(range(len(reflections))).select(z1 == frame - 1)
    shoebox = reflections['shoebox']
    for i in indices2:
      shoebox[i].allocate()
    subset = reflections.select(indices2)
    if len(subset) > 0:
      subset.compute_mask(experiments)

    # Get the mask and data
    raw_mask = imageset.get_mask(frame)[0]
    data = imageset.get_raw_data(frame)[0]

    sbox_mask = reflections.select(indices)['shoebox'].apply_background_mask(frame, 1, (height, width))

    mask = raw_mask & sbox_mask

    from dials.algorithms.image.threshold import DispersionThreshold
    threshold = DispersionThreshold(
      data.all(),
      (3,3),
      6,3,0,2)
    new_mask = flex.bool(mask.accessor())
    threshold(data, mask, new_mask)
    mask = mask & (~new_mask)

    scale_data, scale_mask = fit(data, mask, model_data, model_mask, (9,9))

    fill_mask = scale_mask | ~imageset.get_mask(frame)[0]

    #from matplotlib import pylab
    #pylab.imshow(scale_mask.as_numpy_array(), vmin=0, vmax=2, interpolation='none')
    #pylab.show()
    #pylab.imshow(imageset.get_mask(frame)[0].as_numpy_array(), vmin=0, vmax=2, interpolation='none')
    #pylab.show()
    #pylab.imshow(fill_mask.as_numpy_array(), vmin=0, vmax=2, interpolation='none')
    #pylab.show()

    indices4 = flex.size_t(range(len(scale_data))).select(~fill_mask.as_1d())
    print len(indices4)
    if len(indices4) > 0:
      if indices4[0] == 0:
        indices4 = indices4[1:]
      for i in indices4:
        scale_data[i] = scale_data[i-1]
        scale_mask[i] = scale_mask[i-1]

    #import matplotlib
    #matplotlib.use("Agg")
    #from matplotlib import pylab
    #fig = pylab.figure(dpi=300)
    #pylab.imshow(scale_data.as_numpy_array(), vmin=0, vmax=2, interpolation='none')
    #pylab.colorbar()
    #pylab.savefig("scale_%d.png" % frame)
    #pylab.clf()
    ##pylab.show()
    #exit(0)

    #pylab.hist(scale_data.as_1d().select(scale_mask.as_1d()).as_numpy_array(),
    #           bins=100)
    #pylab.show()

    sd1 = scale_data.as_1d()
    sm1 = scale_mask.as_1d()
    scale_min = flex.min(sd1.select(sm1))
    scale_max = flex.max(sd1.select(sm1))
    scale_avr = flex.sum(sd1.select(sm1)) / sm1.count(True)

    background = model_data * scale_data

    reflections['shoebox'].select(indices).apply_pixel_data(
      data.as_double(),
      background,
      raw_mask,
      frame,
      1)

    subset = reflections.select(indices3)
    if len(subset) > 0:
      subset.compute_summed_intensity()
      subset.compute_centroid(experiments)
      reflections.set_selected(indices3, subset)
    for i in indices3:
      shoebox[i].deallocate()
    num_integrated = subset.get_flags(subset.flags.integrated_sum).count(True)

    print "Image %d: integrated %d reflections" % (frame, num_integrated)
    # print "Image %d: selected %d reflections, scale(min,max_avr)=%f, %f, %f" % (
    #   frame,
    #   len(subset),
    #   scale_min,
    #   scale_max,
    #   scale_avr)


  reflections.as_pickle("reintegrated.pickle")

    #with open("scale_%d.pickle" % frame, "w") as outfile:
    #  pickle.dump((scale_data, scale_mask), outfile, protocol=pickle.HIGHEST_PROTOCOL)
