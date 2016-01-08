

def diffuse(image, mask, max_iter=10, tolerance=1e-5):
  from dials.array_family import flex
  for num_iter in range(max_iter):
    new_image = flex.double(list(image))
    new_image.reshape(image.accessor())
    diff = 0
    num = 0
    print mask.count(1), mask.count(0)
    for j in range(image.all()[0]):
      for i in range(image.all()[1]):
        if mask[j,i] == False:
          if j > 0 and j < image.all()[0]-1 and i > 0 and i < image.all()[1]-1:
            new_image[j,i] = sum([
              image[j-1,i],
              image[j+1,i],
              image[j,i-1],
              image[j,i+1]]) / 4.0
            diff += abs(new_image[j,i] - image[j,i])
            num+=1
    print num_iter, diff / num
    if diff / num < tolerance:
      break
    image = new_image
  return new_image



if __name__ == '__main__':
  import sys
  from dials.algorithms.image.filter import mean_filter
  from dials.array_family import flex

  import pickle

  avr, mask_orig, stddev = pickle.load(open(sys.argv[1]))

  mask = mask_orig.as_1d().as_int()
  #mask = mask * (avr < 750).as_1d().as_int() # 4EPJ
  #mask = mask * (avr < 220).as_1d().as_int() # 4FMR
  mask = mask * (avr < 120).as_1d().as_int() # 4PUC

  indices = flex.size_t(range(len(mask))).select(mask == 0)
  avr2 = avr.as_1d()
  avr2.set_selected(indices, flex.double(len(indices), 0))
  avr2.reshape(avr.accessor())
  avr = avr2
  std2 = stddev.as_1d()
  std2.set_selected(indices, flex.double(len(indices), 0))
  std2.reshape(stddev.accessor())
  stddev = std2
  iod = flex.double(len(avr))
  mask2 = avr2.as_1d() > 0
  indices2 = flex.size_t(range(len(mask))).select(mask2)
  avr2 = avr2.as_1d().select(mask2)
  var2 = std2.as_1d().select(mask2)**2
  iod.set_selected(indices2, var2 / avr2)
  iod.reshape(avr.accessor())
  mask.reshape(avr.accessor())

  print flex.max(avr), flex.min(avr)
  from matplotlib import pylab
  pylab.imshow(iod.as_numpy_array(), vmax=3, interpolation='none')
  pylab.colorbar()
  pylab.show()
  exit(0)
  from matplotlib import pylab
  pylab.imshow(avr.as_numpy_array(), interpolation='none')
  pylab.colorbar()
  pylab.show()
  pylab.hist(avr, bins=100)
  pylab.show()

  #avr = mean_filter(avr, mask, (1,1), 2)
  #print flex.max(avr), flex.min(avr)


  #from matplotlib import pylab
  #pylab.imshow(avr.as_numpy_array(), interpolation='none')
  #pylab.show()

  import cPickle as pickle
  pickle.dump((avr, mask_orig), open("model_filtered.pickle", "w"))

  # avr = avr.as_1d().select(avr.as_1d() > 0.6)

  # pylab.hist(avr, bins=100)
  # pylab.show()

  #from matplotlib import pylab
  #pylab.imshow(mask.as_numpy_array(), vmax=5, interpolation='none')
  #pylab.show()
