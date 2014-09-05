
# def run_test_single(imageset, size):
#   from dials.array_family import flex
#   from time import time
#   import sys
#   print "Allocating: ", size
#   a = flex.double(size)
#   st = time()
#   for image in imageset:
#     sys.stdout.write(".")
#     sys.stdout.flush()
#   print ""

#   ft = time()
#   return size, ft - st

# def run_test(imageset):

#   results = [
#     run_test_single(imageset, 1000000),
#     run_test_single(imageset, 5000000),
#     run_test_single(imageset, 10000000),
#     run_test_single(imageset, 50000000),
#     run_test_single(imageset, 100000000),
#     run_test_single(imageset, 250000000)
#   ]

#   outfile = open("profile.txt", "w")
#   for r in results:
#     print >>outfile, r[0], r[1]


if __name__ == '__main__':
  import sys
  from dxtbx.model.experiment.experiment_list import ExperimentListFactory

  exlist = ExperimentListFactory.from_json_file(sys.argv[1])

  if len(sys.argv) > 2:
    imageset = exlist[0].imageset[0:int(sys.argv[2])]
  else:
    imageset = exlist[0].imageset

  print len(imageset)

  # run_test(imageset)

  from time import time
  st = time()
  n = int(len(imageset) / 4)
  for i in range(0, n):
    image = imageset[i]
    sys.stdout.write('.')
    sys.stdout.flush()
  for i in range(n-17, 2*n):
    image = imageset[i]
    sys.stdout.write('.')
    sys.stdout.flush()
  for i in range(2*n-17, 3*n):
    image = imageset[i]
    sys.stdout.write('.')
    sys.stdout.flush()
  for i in range(3*n-17,len(imageset)):
    image = imageset[i]
    sys.stdout.write('.')
    sys.stdout.flush()
  print ''
  ft = time()
  print ft - st
