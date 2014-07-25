
from __future__ import division
import boost.python

class FastIntegrator(object):

  def __init__(self, exlist, predictions, 
               num_proc=1, mp_method="multiprocessing"):
    from dials.algorithms.integration import FastIntegratorInternal
  
    # Get the imagesets and ensure only 1
    imagesets = exlist.imagesets()
    assert(len(imagesets) == 1)

    # Save some stuff
    self._imageset = imagesets[0]
    self._mp_method = mp_method

    # Get the image range
    image0, image1 = self._imageset.get_scan().get_image_range()

    # Create the integrator
    self._integrator = FastIntegratorInternal(
      predictions, 
      image0,
      image1,
      num_proc)

  def integrate(self):
    from libtbx import easy_mp

    class Worker(object):
      
      def __init__(self, imageset):
        self.imageset = imageset

      def __call__(self, worker):
        from dials.model.data import Image
        index0 = worker.first()
        index1 = worker.last()
        for index in range(index0, index1):
          worker.next(Image(self.imageset[index]))
        return worker.result()

    # Perform the integration on multiple threads
    num_proc = len(self._integrator)
    if num_proc > 1:
      worker_results = easy_mp.parallel_map(
        func=Worker(self._imageset),
        iterable=list(self._integrator.workers()),
        processes=num_proc,
        method=self._mp_method,
        preserve_order=True)
    else:
      worker = Worker(self._imageset)
      worker_results = [worker(self._integrator.worker(0))]

    # Accumulate the results in the worker
    for result in worker_results:
      self._integrator.accumulate(result)

    # Return the result
    return self._integrator.result()
    


def integrate_quickly(exlist, predictions, 
                      num_proc=1, mp_method="multiprocessing"):

  integrator = FastIntegrator(exlist, predictions, num_proc, mp_method)
  return integrator.integrate()
