# LIBTBX_SET_DISPATCHER_NAME dev.dials.scale_down_images

from __future__ import absolute_import, division

def nproc():
  from libtbx.introspection import number_of_processors
  return number_of_processors(return_value_if_unknown=-1)

def joiner(args):
  from dials.util.scale_down_image import scale_down_image
  scale_down_image(*args)
  print args[1]

def scale_down_images(in_template, out_template, start, end, scale_factor):
  from multiprocessing import Pool

  jobs = [(in_template % j, out_template % j, scale_factor) for j in
          range(start, end + 1)]

  pool = Pool(processes=nproc())
  result = pool.map_async(joiner, jobs)
  result.get()

  return result

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 6:
    raise RuntimeError, '%s in_\%04d.cbf out_\%04d.cbf start end scale' % \
      sys.argv[0]
  in_template = sys.argv[1]
  out_template = sys.argv[2]
  start = int(sys.argv[3])
  end = int(sys.argv[4])
  scale_factor = float(sys.argv[5])
  scale_down_images(in_template, out_template, start, end, scale_factor)
