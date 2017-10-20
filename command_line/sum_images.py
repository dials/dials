# LIBTBX_SET_DISPATCHER_NAME dev.dials.sum_images

from __future__ import absolute_import, division

def sum_images(in_template, out_image, start, end):
  from dials.util.rebin_images import main_sum
  in_images = [in_template % j for j in range(start, end + 1)]
  main_sum(in_images, out_image)
  return

if __name__ == '__main__':
  import sys
  if len(sys.argv) != 5:
    raise RuntimeError('%s in_\%d_0001.cbf out_1_0001.cbf start end' % \
      sys.argv[0])
  in_template = sys.argv[1]
  out_image = sys.argv[2]
  start = int(sys.argv[3])
  end = int(sys.argv[4])
  sum_images(in_template, out_image, start, end)
