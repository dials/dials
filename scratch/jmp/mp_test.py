


def test_function2():
  from dials.util.options import OptionParser
  from time import time
  from omptbx import omp_get_max_threads


  args = [
    '/media/upc86896/Betelgeuse/Data/xia2_test_data/i19_sucrose/processed/dataset1/experiments.json',
    '/media/upc86896/Betelgeuse/Data/xia2_test_data/i19_sucrose/processed/dataset1/shoeboxes_0.pickle'
  ]

  parser = OptionParser(read_experiments=True, read_reflections=True)
  params, options =  parser.parse_args(args=args)

  experiment = params.input.experiments[0].data[0]
  reflections = params.input.reflections[0].data

  print "N Threads: ", omp_get_max_threads()
  print "N Refl: ", len(reflections)

  from dials.algorithms.integration.fitrs import test_function

  st = time()
  test_function(
    experiment.beam,
    experiment.detector,
    experiment.goniometer,
    experiment.scan,
    0.071016,
    0.390601,
    5,
    reflections)
  print "Time: ", time() - st

if __name__ == '__main__':

  test_function2()
